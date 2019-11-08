import os
from ROOT import RooFit as rf
from ROOT import TFile, TTree, TCanvas
from ROOT import RooWorkspace, RooDataSet, RooArgSet, RooSimultaneous, RooCmdArg, RooKeysPdf
from utils import *

path = "/Users/matt/Scratch/lhcb/ArnauFits/Tuples/"

### SIGNAL SHAPES ###
# These are fully reco'd and should be fairly similar shapes
# Require them to have the same tail params and CB fraction
# Allow a shift in mean and a scale in width for g -> pi0
# Allow a shift in mean and a scale in width for b0 -> bs0
def SignalShape(free=True,sim=True):

  w = RooWorkspace('w_sig','w_sig')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))


  files = { 'b0g'   : path+"B0gamma/9_MC_B0gamma_Signal_BestCut.root",
            'b0pi0' : path+"B0pi0/9_MC_B0pi0_Signal_BestCut.root",
            'bsg'   : path+"Bsgamma/9_MC_Bsgamma_Signal_BestCut.root",
            'bspi0' : path+"Bspi0/9_MC_Bspi0_Signal_BestCut.root"
          }

  ### Make the MC dsets ###
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("sig_mc_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  ### Make the signal PDF ###
  SharedParams(w)

  # B0g signal model
  w.factory("b0g_mean[5279.63,5250,5400]")
  w.factory("b0g_sigma[30,0,50]")
  w.factory("alpha1[2,0,4]")
  w.factory("n1[3]")
  w.factory("alpha2[-2,-4,0]")
  w.factory("n2[3]")
  w.factory("f1[0.5,0.,1.]")
  w.factory("CBShape::b0g_cb1( B_DTFDict_D0_B_M, b0g_mean, b0g_sigma, alpha1, n1 )")
  w.factory("CBShape::b0g_cb2( B_DTFDict_D0_B_M, b0g_mean, b0g_sigma, alpha2, n2 )")
  w.factory("SUM::sig_mc_pdf_b0g( f1*b0g_cb1, b0g_cb2 )")
  # create dm shift and sigma scale for g->pi0
  w.factory("sum::b0pi0_mean(b0g_mean,dm_g2pi0)")
  w.factory("prod::b0pi0_sigma(b0g_sigma,ssig_g2pi0)")
  # B0pi0 signal model
  w.factory("CBShape::b0pi0_cb1( B_DTFDict_D0_B_M, b0pi0_mean, b0pi0_sigma, alpha1, n1 )")
  w.factory("CBShape::b0pi0_cb2( B_DTFDict_D0_B_M, b0pi0_mean, b0pi0_sigma, alpha2, n2 )")
  w.factory("SUM::sig_mc_pdf_b0pi0( f1*b0pi0_cb1, b0pi0_cb2 )")
  # create dm shift and sigma scale for b0 -> bs
  w.factory("sum::bsg_mean(b0g_mean,dm_b02bs)")
  w.factory("prod::bsg_sigma(b0g_sigma,ssig_b02bs)")
  # Bsg signal model
  w.factory("CBShape::bsg_cb1( B_DTFDict_D0_B_M, bsg_mean, bsg_sigma, alpha1, n1 )")
  w.factory("CBShape::bsg_cb2( B_DTFDict_D0_B_M, bsg_mean, bsg_sigma, alpha2, n2 )")
  w.factory("SUM::sig_mc_pdf_bsg( f1*bsg_cb1, bsg_cb2 )")
  # create dm shift and sigma scale for b0g -> bspi0
  w.factory("sum::bspi0_mean(b0g_mean,dm_g2pi0,dm_b02bs)")
  w.factory("prod::bspi0_sigma(b0g_sigma,ssig_g2pi0,ssig_b02bs)")
  # Bspi0 signal model
  w.factory("CBShape::bspi0_cb1( B_DTFDict_D0_B_M, bspi0_mean, bspi0_sigma, alpha1, n1 )")
  w.factory("CBShape::bspi0_cb2( B_DTFDict_D0_B_M, bspi0_mean, bspi0_sigma, alpha2, n2 )")
  w.factory("SUM::sig_mc_pdf_bspi0( f1*bspi0_cb1, bspi0_cb2 )")

  # make sim pdf and sim data
  CreateSimPdf(w,'sig_mc')
  CreateSimData(w,'sig_mc')

  c = TCanvas('c_sig_mc','c',1400,1200)
  c.Divide(2,2)

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'sig_mc_pdf_%s'%(samp)
      dsetname = 'sig_mc_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname))
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp)
    c.Update()
    c.Modified()
    c.Print("plots/sig_mc_free_fit.pdf")

  if sim:
    pdfname = 'sig_mc_sim_pdf'
    dsetname = 'sig_mc_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    for i, samp in enumerate(samples):
      dopts = RooLinkedList()
      dopts.Add( RooCmdArg( rf.Cut('cat==cat::%s'%samp) ) )
      popts = RooLinkedList()
      popts.Add( RooCmdArg( rf.Slice(w.cat("cat"),samp) ) )
      if w.data(dsetname): popts.Add( RooCmdArg( rf.ProjWData(w.data(dsetname)) ) )
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp, popts, dopts)

    c.Update()
    c.Modified()
    c.Print("plots/sig_mc_sim_fit.pdf")

    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.saveSnapshot(pdfname.replace('pdf','fit'), pars)

  w.writeToFile('files/w_sig.root')

### Mis Reco shape ###
### These will be similar for B0 -> Bs (with a shift in mean)
### But rather different for g -> pi0 because in the first case it's pi0 with a g missing
### in the second case it's a g with an extra fake or real g added
def MisRecShape(free=True,sim=True):

  w = RooWorkspace('w_misrec','w_misrec')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))

  files = { 'b0g'   : path+"B0gamma/9_MC_B0gamma_MisSignal_BestCut.root",
            'b0pi0' : path+"B0pi0/9_MC_B0pi0_AlmostSignal_pi0s_BestCut.root",
            'bsg'   : path+"Bsgamma/9_MC_Bsgamma_MisSignal_Matching_BestCut.root",
            'bspi0' : path+"Bspi0/9_MC_Bspi0_AlmostSignal_BestCut.root"
          }

  ### Make the MC dsets ###
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("misrec_mc_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  ### Make the misrec PDF ###
  SharedParams(w)

  ## parameters
  tfs = TFile('files/w_sig.root')
  ws = tfs.Get('w_sig')
  getattr(w,'import')(ws.var('b0g_mean'))
  getattr(w,'import')(ws.function('b0pi0_mean'))
  getattr(w,'import')(ws.function('bsg_mean'))
  getattr(w,'import')(ws.function('bspi0_mean'))
  tfs.Close()
  w.var('b0g_mean').setConstant(True)
  w.var('dm_b02bs').setConstant(True)
  w.var('dm_g2pi0').setConstant(True)

  # means (should be relative to signal mean)
  w.factory("dm_misrec2sig_b0g[10,-50,50]")
  w.factory("dm_misrec2sig_b0pi0[10,-50,50]")
  w.factory("dm_misrec2sig_bsg[10,-50,50]")
  w.factory("dm_misrec2sig_bspi0[10,-50,50]")
  w.factory("sum::b0g_misrec_mean( dm_misrec2sig_b0g, b0g_mean )")
  w.factory("sum::b0pi0_misrec_mean( dm_misrec2sig_b0pi0, b0pi0_mean )")
  w.factory("sum::bsg_misrec_mean( dm_misrec2sig_bsg, bsg_mean )")
  w.factory("sum::bspi0_misrec_mean( dm_misrec2sig_bspi0, bspi0_mean )")

  # sigmas
  w.factory("b0g_misrec_sigma[30,0,100]")
  w.factory("prod::bsg_misrec_sigma(b0g_misrec_sigma,ssig_b02bs)")
  # a scale relating pi0 with g to g with extra g
  w.factory("ssig_missg2addg[1,0.2,2]")
  w.factory("prod::b0pi0_misrec_sigma( b0g_misrec_sigma, ssig_missg2addg )")
  w.factory("prod::bspi0_misrec_sigma( b0g_misrec_sigma, ssig_b02bs, ssig_missg2addg )")
  w.factory("b0pi0_misrec_sigma2[80,0,200]")
  w.factory("prod::bspi0_misrec_sigma2( b0pi0_misrec_sigma2, ssig_b02bs )")

  # tails params
  w.factory("misrec_alpha1_g[2.1,0,4]")
  w.factory("misrec_alpha1_pi0[2.1,0,4]")
  w.factory("misrec_n1[3]")
  w.factory("misrec_alpha2_g[-2.1,-4,0]")
  w.factory("misrec_alpha2_pi0[-2.1,-4,0]")
  w.factory("misrec_n2[3]")
  w.factory("misrec_f1_g[0.5,0.,1]")
  w.factory("misrec_f1_pi0[0.5,0.,1]")

  # shapes
  # B0g model
  w.factory("CBShape::misrec_b0g_cb1( B_DTFDict_D0_B_M, b0g_misrec_mean, b0g_misrec_sigma, misrec_alpha1_g, misrec_n2)")
  w.factory("CBShape::misrec_b0g_cb2( B_DTFDict_D0_B_M, b0g_misrec_mean, b0g_misrec_sigma, misrec_alpha2_g, misrec_n2)")
  w.factory("SUM::misrec_mc_pdf_b0g( misrec_f1_g*misrec_b0g_cb1, misrec_b0g_cb2 )")
  # B0pi0  model
  w.factory("CBShape::misrec_b0pi0_cb1( B_DTFDict_D0_B_M, b0pi0_misrec_mean, b0pi0_misrec_sigma, misrec_alpha1_pi0, misrec_n1)")
  w.factory("CBShape::misrec_b0pi0_cb2( B_DTFDict_D0_B_M, b0pi0_misrec_mean, b0pi0_misrec_sigma2, misrec_alpha2_pi0, misrec_n2)")
  w.factory("SUM::misrec_mc_pdf_b0pi0( misrec_f1_pi0*misrec_b0pi0_cb1, misrec_b0pi0_cb2 )")
  # Bsg  model
  w.factory("CBShape::misrec_bsg_cb1( B_DTFDict_D0_B_M, bsg_misrec_mean, bsg_misrec_sigma, misrec_alpha1_g, misrec_n1)")
  w.factory("CBShape::misrec_bsg_cb2( B_DTFDict_D0_B_M, bsg_misrec_mean, bsg_misrec_sigma, misrec_alpha2_g, misrec_n2)")
  w.factory("SUM::misrec_mc_pdf_bsg( misrec_f1_g*misrec_bsg_cb1, misrec_bsg_cb2 )")
  # Bspi0  model
  w.factory("CBShape::misrec_bspi0_cb1( B_DTFDict_D0_B_M, bspi0_misrec_mean, bspi0_misrec_sigma, misrec_alpha1_pi0, misrec_n1)")
  w.factory("CBShape::misrec_bspi0_cb2( B_DTFDict_D0_B_M, bspi0_misrec_mean, bspi0_misrec_sigma2, misrec_alpha2_pi0, misrec_n2)")
  w.factory("SUM::misrec_mc_pdf_bspi0( misrec_f1_pi0*misrec_bspi0_cb1, misrec_bspi0_cb2 )")

  # make sim pdf and sim data
  CreateSimPdf(w,'misrec_mc')
  CreateSimData(w,'misrec_mc')

  c = TCanvas('c_sig_mc','c',1400,1200)
  c.Divide(2,2)

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'misrec_mc_pdf_%s'%(samp)
      dsetname = 'misrec_mc_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname))
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp)
    c.Update()
    c.Modified()
    c.Print("plots/misrec_mc_free_fit.pdf")

  if sim:
    pdfname = 'misrec_mc_sim_pdf'
    dsetname = 'misrec_mc_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    for i, samp in enumerate(samples):
      dopts = RooLinkedList()
      dopts.Add( RooCmdArg( rf.Cut('cat==cat::%s'%samp) ) )
      popts = RooLinkedList()
      popts.Add( RooCmdArg( rf.Slice(w.cat("cat"),samp) ) )
      if w.data(dsetname): popts.Add( RooCmdArg( rf.ProjWData(w.data(dsetname)) ) )
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp, popts, dopts)

    c.Update()
    c.Modified()
    c.Print("plots/misrec_mc_sim_fit.pdf")

    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.saveSnapshot(pdfname.replace('pdf','fit'), pars)

  w.writeToFile('files/w_misrec.root')

### B -> D*pipi shape ###
### should be identical for B0 and Bs because for both you have a K misID's as pi
### differences in terms of shifts and scales for g -> pi0 should be the same as for
### signal and other parts
def BDstPPShape(free=True,sim=True):
  w = RooWorkspace('w_bdstpp','w_bdstpp')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))

  files = { 'b0g'   : path+"B0gamma/9_MC_B0gamma_Bg_B2Dstpipi_BestCut.root",
            'b0pi0' : path+"B0pi0/9_MC_B0pi0_Bg_B2Dstpipi_BestCut.root",
            'bsg'   : path+"Bsgamma/9_MC_Bsgamma_Bg_B2Dstpipi_BestCut.root",
            'bspi0' : path+"Bspi0/9_MC_Bspi0_Bg_B2Dstpipi_BestCut.root"
          }

  ### Make the MC dsets ###
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("bdstpp_mc_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  ### Make the signal PDF ###
  SharedParams(w)

  # parameters
  tfs = TFile('files/w_sig.root')
  ws = tfs.Get('w_sig')
  getattr(w,'import')(ws.var('b0g_mean'))
  getattr(w,'import')(ws.function('b0pi0_mean'))
  getattr(w,'import')(ws.function('bsg_mean'))
  getattr(w,'import')(ws.function('bspi0_mean'))
  tfs.Close()
  w.var('b0g_mean').setConstant(True)
  w.var('dm_b02bs').setConstant(True)
  w.var('dm_g2pi0').setConstant(True)

  # means (should be relative to signal mean)
  w.factory("dm_bdstpp2sig_b0g[10,-50,50]")
  w.factory("dm_bdstpp2sig_b0pi0[10,-50,50]")
  w.factory("sum::bdstpp_b0g_mean( dm_bdstpp2sig_b0g, b0g_mean )")
  w.factory("sum::bdstpp_b0pi0_mean( dm_bdstpp2sig_b0pi0, b0pi0_mean )")

  # sigmas
  w.factory("bdstpp_b0g_sigma[30,0,50]")
  w.factory("prod::bdstpp_b0pi0_sigma(bdstpp_b0g_sigma,ssig_g2pi0)")

  # tail params
  w.factory("bdstpp_alpha1[2,0,4]")
  w.factory("bdstpp_n1[3]")
  w.factory("bdstpp_alpha2[-2,-4,0]")
  w.factory("bdstpp_n2[3]")
  w.factory("bdstpp_f1[0.5,0.,1.]")

  # B0g signal model
  w.factory("CBShape::bdstpp_b0g_cb1( B_DTFDict_D0_B_M, bdstpp_b0g_mean, bdstpp_b0g_sigma, bdstpp_alpha1, bdstpp_n1 )")
  w.factory("CBShape::bdstpp_b0g_cb2( B_DTFDict_D0_B_M, bdstpp_b0g_mean, bdstpp_b0g_sigma, bdstpp_alpha2, bdstpp_n2 )")
  w.factory("SUM::bdstpp_mc_pdf_b0g( bdstpp_f1*bdstpp_b0g_cb1, bdstpp_b0g_cb2 )")
  # B0pi0 signal model
  w.factory("CBShape::bdstpp_b0pi0_cb1( B_DTFDict_D0_B_M, bdstpp_b0pi0_mean, bdstpp_b0pi0_sigma, bdstpp_alpha1, bdstpp_n1 )")
  w.factory("CBShape::bdstpp_b0pi0_cb2( B_DTFDict_D0_B_M, bdstpp_b0pi0_mean, bdstpp_b0pi0_sigma, bdstpp_alpha2, bdstpp_n2 )")
  w.factory("SUM::bdstpp_mc_pdf_b0pi0( bdstpp_f1*bdstpp_b0pi0_cb1, bdstpp_b0pi0_cb2 )")
  # Bsg signal model
  w.factory("CBShape::bdstpp_bsg_cb1( B_DTFDict_D0_B_M, bdstpp_b0g_mean, bdstpp_b0g_sigma, bdstpp_alpha1, bdstpp_n1 )")
  w.factory("CBShape::bdstpp_bsg_cb2( B_DTFDict_D0_B_M, bdstpp_b0g_mean, bdstpp_b0g_sigma, bdstpp_alpha2, bdstpp_n2 )")
  w.factory("SUM::bdstpp_mc_pdf_bsg( bdstpp_f1*bdstpp_bsg_cb1, bdstpp_bsg_cb2 )")
  # Bspi0 signal model
  w.factory("CBShape::bdstpp_bspi0_cb1( B_DTFDict_D0_B_M, bdstpp_b0pi0_mean, bdstpp_b0pi0_sigma, bdstpp_alpha1, bdstpp_n1 )")
  w.factory("CBShape::bdstpp_bspi0_cb2( B_DTFDict_D0_B_M, bdstpp_b0pi0_mean, bdstpp_b0pi0_sigma, bdstpp_alpha2, bdstpp_n2 )")
  w.factory("SUM::bdstpp_mc_pdf_bspi0( bdstpp_f1*bdstpp_bspi0_cb1, bdstpp_bspi0_cb2 )")

  # make sim pdf and sim data
  CreateSimPdf(w,'bdstpp_mc')
  CreateSimData(w,'bdstpp_mc')

  c = TCanvas('c_bdstpp_mc','c',1400,1200)
  c.Divide(2,2)

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'bdstpp_mc_pdf_%s'%(samp)
      dsetname = 'bdstpp_mc_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname))
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp)
    c.Update()
    c.Modified()
    c.Print("plots/bdstpp_mc_free_fit.pdf")

  if sim:
    pdfname = 'bdstpp_mc_sim_pdf'
    dsetname = 'bdstpp_mc_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    for i, samp in enumerate(samples):
      dopts = RooLinkedList()
      dopts.Add( RooCmdArg( rf.Cut('cat==cat::%s'%samp) ) )
      popts = RooLinkedList()
      popts.Add( RooCmdArg( rf.Slice(w.cat("cat"),samp) ) )
      if w.data(dsetname): popts.Add( RooCmdArg( rf.ProjWData(w.data(dsetname)) ) )
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp, popts, dopts)

    c.Update()
    c.Modified()
    c.Print("plots/bdstpp_mc_sim_fit.pdf")

    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.saveSnapshot(pdfname.replace('pdf','fit'), pars)

  w.writeToFile('files/w_bdstpp.root')

### B -> D*KK shape ###
### Have two contributions here from B0 -> D*0 K+K- and Bs0 -> D*0 K+K-
### These should be identical for what we reconstruct as B0 and Bs because
### for both you have a pi misID'd as a K
### the difference between the two samples shift and scales for b0 -> bs0 should be
### the same as for signal and other parts
### differences in terms of shifts and scales for g -> pi0 should be the same as for
### signal and other parts
### Note for this bit we will now have 8 categories to control the shape
### I have suspicion this sample is not truth matched properly (what is the guff below the peaks)?
def BDstKKShape(free=True,sim=True):
  w = RooWorkspace('w_bdstkk','w_bdstkk')

  samples = [ 'b0_b0g', 'b0_b0pi0', 'b0_bsg', 'b0_bspi0',
              'bs_b0g', 'bs_b0pi0', 'bs_bsg', 'bs_bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))

  files = { 'b0_b0g'   : path+"B0gamma/9_MC_B0gamma_Bg_B2DstKK_BestCut.root",
            'b0_b0pi0' : path+"B0pi0/9_MC_B0pi0_Bg_B2DstKK_BestCut.root",
            'b0_bsg'   : path+"Bsgamma/9_MC_Bsgamma_Bg_B2DstKK_BestCut.root",
            'b0_bspi0' : path+"Bspi0/9_MC_Bspi0_Bg_B2DstKK_BestCut.root",
            'bs_b0g'   : path+"B0gamma/9_MC_B0gamma_Bg_Bs2DstKK_BestCut.root",
            'bs_b0pi0' : path+"B0pi0/9_MC_B0pi0_Bg_Bs2DstKK_BstCut.root",
            'bs_bsg'   : path+"Bsgamma/9_MC_Bsgamma_Bg_Bs2DstKK_BestCut.root",
            'bs_bspi0' : path+"Bspi0/9_MC_Bspi0_Bg_Bs2DstKK_BestCut.root"
          }

  ### Make the MC dsets ###
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("bdstkk_mc_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  ### Make the signal PDF ###
  SharedParams(w)

  ## FOR THE B -> D* KK sample
  # B0g signal model
  w.factory("bdstkk_b0_b0g_mean[5279.63,5200,5300]")
  w.factory("bdstkk_b0_b0g_sigma[30,10,50]")
  w.factory("bdstkk_alpha1[2,0,4]")
  w.factory("bdstkk_n1[3]")
  w.factory("bdstkk_alpha2[-2,-4,0]")
  w.factory("bdstkk_n2[3]")
  w.factory("bdstkk_f1[0.5,0.,1.]")
  w.factory("CBShape::bdstkk_b0_b0g_cb1( B_DTFDict_D0_B_M, bdstkk_b0_b0g_mean, bdstkk_b0_b0g_sigma, bdstkk_alpha1, bdstkk_n1 )")
  w.factory("CBShape::bdstkk_b0_b0g_cb2( B_DTFDict_D0_B_M, bdstkk_b0_b0g_mean, bdstkk_b0_b0g_sigma, bdstkk_alpha2, bdstkk_n2 )")
  w.factory("SUM::bdstkk_b0_b0g_cbs( bdstkk_f1*bdstkk_b0_b0g_cb1, bdstkk_b0_b0g_cb2 )" )
  #w.factory("Gaussian::bdstkk_b0_b0g_gaus( B_DTFDict_D0_B_M, bdstkk_b0_b0g_mean, bdstkk_b0_b0g_sigma )")
  w.factory("bdstkk_b0_m0[5170,5140,5280]")
  w.factory("bdstkk_b0_c[1,0,40]")
  w.factory("ArgusBG::bdstkk_b0_b0g_bg( B_DTFDict_D0_B_M, bdstkk_b0_m0, bdstkk_b0_c )")
  w.factory("bdstkk_b0_b0g_f2[0.5,0.5,1.]")
  w.factory("SUM::bdstkk_mc_pdf_b0_b0g( bdstkk_b0_b0g_f2*bdstkk_b0_b0g_cbs, bdstkk_b0_b0g_bg )")
  # create dm shift and sigma scale for g->pi0
  w.factory("sum::bdstkk_b0_b0pi0_mean(bdstkk_b0_b0g_mean,dm_g2pi0)")
  w.factory("prod::bdstkk_b0_b0pi0_sigma(bdstkk_b0_b0g_sigma,ssig_g2pi0)")
  # B0pi0 signal model
  w.factory("CBShape::bdstkk_b0_b0pi0_cb1( B_DTFDict_D0_B_M, bdstkk_b0_b0pi0_mean, bdstkk_b0_b0pi0_sigma, bdstkk_alpha1, bdstkk_n1 )")
  w.factory("CBShape::bdstkk_b0_b0pi0_cb2( B_DTFDict_D0_B_M, bdstkk_b0_b0pi0_mean, bdstkk_b0_b0pi0_sigma, bdstkk_alpha2, bdstkk_n2 )")
  w.factory("SUM::bdstkk_b0_b0pi0_cbs( bdstkk_f1*bdstkk_b0_b0pi0_cb1, bdstkk_b0_b0pi0_cb2 )")
  #w.factory("Gaussian::bdstkk_b0_b0pi0_gaus( B_DTFDict_D0_B_M, bdstkk_b0_b0pi0_mean, bdstkk_b0_b0pi0_sigma )")
  w.factory("bdstkk_b0_b0pi0_f2[0.5,0.,1.]")
  w.factory("SUM::bdstkk_mc_pdf_b0_b0pi0( bdstkk_b0_b0pi0_f2*bdstkk_b0_b0pi0_cbs, bdstkk_b0_b0g_bg )")
  # Bsg signal model
  w.factory("CBShape::bdstkk_b0_bsg_cb1( B_DTFDict_D0_B_M, bdstkk_b0_b0g_mean, bdstkk_b0_b0g_sigma, bdstkk_alpha1, bdstkk_n1 )")
  w.factory("CBShape::bdstkk_b0_bsg_cb2( B_DTFDict_D0_B_M, bdstkk_b0_b0g_mean, bdstkk_b0_b0g_sigma, bdstkk_alpha2, bdstkk_n2 )")
  w.factory("SUM::bdstkk_b0_bsg_cbs( bdstkk_f1*bdstkk_b0_bsg_cb1, bdstkk_b0_bsg_cb2 )")
  #w.factory("Gaussian::bdstkk_b0_bsg_gaus( B_DTFDict_D0_B_M, bdstkk_b0_b0g_mean, bdstkk_b0_b0g_sigma )")
  w.factory("SUM::bdstkk_mc_pdf_b0_bsg( bdstkk_b0_b0g_f2*bdstkk_b0_bsg_cbs, bdstkk_b0_b0g_bg )")
  # Bspi0 signal model
  w.factory("CBShape::bdstkk_b0_bspi0_cb1( B_DTFDict_D0_B_M, bdstkk_b0_b0pi0_mean, bdstkk_b0_b0pi0_sigma, bdstkk_alpha1, bdstkk_n1 )")
  w.factory("CBShape::bdstkk_b0_bspi0_cb2( B_DTFDict_D0_B_M, bdstkk_b0_b0pi0_mean, bdstkk_b0_b0pi0_sigma, bdstkk_alpha2, bdstkk_n2 )")
  w.factory("SUM::bdstkk_b0_bspi0_cbs( bdstkk_f1*bdstkk_b0_bspi0_cb1, bdstkk_b0_bspi0_cb2 )")
  #w.factory("Gaussian::bdstkk_b0_bspi0_gaus( B_DTFDict_D0_B_M, bdstkk_b0_b0pi0_mean, bdstkk_b0_b0pi0_sigma )")
  w.factory("SUM::bdstkk_mc_pdf_b0_bspi0( bdstkk_b0_b0pi0_f2*bdstkk_b0_bspi0_cbs, bdstkk_b0_b0g_bg )")

  ## FOR THE Bs -> D* KK sample
  # B0g signal model
  w.factory("sum::bdstkk_bs_b0g_mean( bdstkk_b0_b0g_mean, dm_b02bs )")
  w.factory("prod::bdstkk_bs_b0g_sigma( bdstkk_b0_b0g_sigma, ssig_b02bs )")
  #w.factory("bdstkk_alpha1[2,0,4]")
  #w.factory("bdstkk_n1[3]")
  #w.factory("bdstkk_alpha2[-2,-4,0]")
  #w.factory("bdstkk_n2[3]")
  #w.factory("bdstkk_f1[0.5,0.,1.]")
  w.factory("CBShape::bdstkk_bs_b0g_cb1( B_DTFDict_D0_B_M, bdstkk_bs_b0g_mean, bdstkk_bs_b0g_sigma, bdstkk_alpha1, bdstkk_n1 )")
  w.factory("CBShape::bdstkk_bs_b0g_cb2( B_DTFDict_D0_B_M, bdstkk_bs_b0g_mean, bdstkk_bs_b0g_sigma, bdstkk_alpha2, bdstkk_n2 )")
  w.factory("SUM::bdstkk_bs_b0g_cbs( bdstkk_f1*bdstkk_bs_b0g_cb1, bdstkk_bs_b0g_cb2 )")
  #w.factory("Gaussian::bdstkk_bs_b0g_gaus( B_DTFDict_D0_B_M, bdstkk_bs_b0g_mean, bdstkk_bs_b0g_sigma )")
  #w.factory("bdstkk_bs_m0[5230,5220,5350]")
  w.factory("sum::bdstkk_bs_m0( bdstkk_b0_m0, dm_b02bs )")
  w.factory("bdstkk_bs_c[1,0.01,20]")
  w.factory("ArgusBG::bdstkk_bs_b0g_bg( B_DTFDict_D0_B_M, bdstkk_bs_m0, bdstkk_bs_c )")
  w.factory("bdstkk_bs_b0g_f2[0.5,0.,1.]")
  w.factory("SUM::bdstkk_mc_pdf_bs_b0g( bdstkk_bs_b0g_f2*bdstkk_bs_b0g_cbs, bdstkk_bs_b0g_bg )")
  # create dm shift and sigma scale for g->pi0
  w.factory("sum::bdstkk_bs_b0pi0_mean(bdstkk_bs_b0g_mean,dm_g2pi0)")
  w.factory("prod::bdstkk_bs_b0pi0_sigma(bdstkk_bs_b0g_sigma,ssig_g2pi0)")
  # B0pi0 signal model
  w.factory("CBShape::bdstkk_bs_b0pi0_cb1( B_DTFDict_D0_B_M, bdstkk_bs_b0pi0_mean, bdstkk_bs_b0pi0_sigma, bdstkk_alpha1, bdstkk_n1 )")
  w.factory("CBShape::bdstkk_bs_b0pi0_cb2( B_DTFDict_D0_B_M, bdstkk_bs_b0pi0_mean, bdstkk_bs_b0pi0_sigma, bdstkk_alpha2, bdstkk_n2 )")
  w.factory("SUM::bdstkk_bs_b0pi0_cbs( bdstkk_f1*bdstkk_bs_b0pi0_cb1, bdstkk_bs_b0pi0_cb2 )")
  #w.factory("Gaussian::bdstkk_bs_b0pi0_gaus( B_DTFDict_D0_B_M, bdstkk_bs_b0pi0_mean, bdstkk_bs_b0pi0_sigma )")
  w.factory("bdstkk_bs_b0pi0_f2[0.5,0.,1.]")
  w.factory("SUM::bdstkk_mc_pdf_bs_b0pi0( bdstkk_bs_b0pi0_f2*bdstkk_bs_b0pi0_cbs, bdstkk_bs_b0g_bg )")
  # Bsg signal model
  w.factory("CBShape::bdstkk_bs_bsg_cb1( B_DTFDict_D0_B_M, bdstkk_bs_b0g_mean, bdstkk_bs_b0g_sigma, bdstkk_alpha1, bdstkk_n1 )")
  w.factory("CBShape::bdstkk_bs_bsg_cb2( B_DTFDict_D0_B_M, bdstkk_bs_b0g_mean, bdstkk_bs_b0g_sigma, bdstkk_alpha2, bdstkk_n2 )")
  w.factory("SUM::bdstkk_bs_bsg_cbs( bdstkk_f1*bdstkk_bs_bsg_cb1, bdstkk_bs_bsg_cb2 )")
  #w.factory("Gaussian::bdstkk_bs_bsg_gaus( B_DTFDict_D0_B_M, bdstkk_bs_b0g_mean, bdstkk_bs_b0g_sigma )")
  w.factory("SUM::bdstkk_mc_pdf_bs_bsg( bdstkk_bs_b0g_f2*bdstkk_bs_bsg_cbs, bdstkk_bs_b0g_bg )")
  # Bspi0 signal model
  w.factory("CBShape::bdstkk_bs_bspi0_cb1( B_DTFDict_D0_B_M, bdstkk_bs_b0pi0_mean, bdstkk_bs_b0pi0_sigma, bdstkk_alpha1, bdstkk_n1 )")
  w.factory("CBShape::bdstkk_bs_bspi0_cb2( B_DTFDict_D0_B_M, bdstkk_bs_b0pi0_mean, bdstkk_bs_b0pi0_sigma, bdstkk_alpha2, bdstkk_n2 )")
  w.factory("SUM::bdstkk_bs_bspi0_cbs( bdstkk_f1*bdstkk_bs_bspi0_cb1, bdstkk_bs_bspi0_cb2 )")
  #w.factory("Gaussian::bdstkk_bs_bspi0_gaus( B_DTFDict_D0_B_M, bdstkk_bs_b0pi0_mean, bdstkk_bs_b0pi0_sigma )")
  w.factory("SUM::bdstkk_mc_pdf_bs_bspi0( bdstkk_bs_b0pi0_f2*bdstkk_bs_bspi0_cbs, bdstkk_bs_b0g_bg )")

  # make sim pdf and sim data
  CreateSimPdf(w,'bdstkk_mc')
  CreateSimData(w,'bdstkk_mc')

  c1 = TCanvas('c1_bdstkk_mc','c',1400,1200)
  c1.Divide(2,2)
  c2 = TCanvas('c2_bdstkk_mc','c',1400,1200)
  c2.Divide(2,2)

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'bdstkk_mc_pdf_%s'%(samp)
      dsetname = 'bdstkk_mc_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname))
      if i<=3: plot( c1.cd(i+1), w, pdfname, dsetname, True, samp)
      else:    plot( c2.cd(i+1-4), w, pdfname, dsetname, True, samp)
    c1.Update()
    c1.Modified()
    c1.Print("plots/bdstkk_mc_free_fit1.pdf")
    c2.Update()
    c2.Modified()
    c2.Print("plots/bdstkk_mc_free_fit2.pdf")

  if sim:
    pdfname = 'bdstkk_mc_sim_pdf'
    dsetname = 'bdstkk_mc_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    for i, samp in enumerate(samples):
      dopts = RooLinkedList()
      dopts.Add( RooCmdArg( rf.Cut('cat==cat::%s'%samp) ) )
      popts = RooLinkedList()
      popts.Add( RooCmdArg( rf.Slice(w.cat("cat"),samp) ) )
      if w.data(dsetname): popts.Add( RooCmdArg( rf.ProjWData(w.data(dsetname)) ) )
      if i<=3: plot( c1.cd(i+1), w, pdfname, dsetname, True, samp, popts, dopts)
      else:    plot( c2.cd(i+1-4), w, pdfname, dsetname, True, samp, popts, dopts)

    c1.Update()
    c1.Modified()
    c1.Print("plots/bdstkk_mc_sim_fit1.pdf")
    c2.Update()
    c2.Modified()
    c2.Print("plots/bdstkk_mc_sim_fit2.pdf")

    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.saveSnapshot(pdfname.replace('pdf','fit'), pars)

  w.writeToFile('files/w_bdstkk.root')

### B -> DKpi shape ###
### Will be B0 or Bs -> DKpi with random added g or pi0
### So B0 to Bs should have the normal shift and smear
### For g -> pi0 mode the shift and
### smear are probably related to the misrec case but
### we'll treat them as indepdent
### For some reason there's seems to be two components, a small
### one at the right mass and a large one at the wrong mass
### Maybe in one case some FSR or v. soft netutral is added?

def BDKPShape(free=True,sim=True):
  w = RooWorkspace('w_bdkp','w_bdkp')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))

  files = { 'b0g'   : path+"B0gamma/9_MC_B0gamma_Bg_B2DKpi_BestCut.root",
            'b0pi0' : path+"B0pi0/9_MC_B0pi0_Bg_B2DKpi_BestCut.root",
            'bsg'   : path+"Bsgamma/9_MC_Bsgamma_Bg_Bs2DKpi_BestCut.root",
            'bspi0' : path+"Bspi0/9_MC_Bspi0_Bg_Bs2DKpi_BestCut.root"
          }

  ### Make the MC dsets ###
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("bdkp_mc_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  ### Make the signal PDF ###
  SharedParams(w)

  # means
  w.factory("bdkp_b0g_mean[5450,5400,5550]")
  w.factory("sum::bdkp_bsg_mean(bdkp_b0g_mean, dm_b02bs)")
  # dm shift for g -> pi0 in this case
  w.factory("dm_addedg2pi0[-10,-50,50]")
  w.factory("sum::bdkp_b0pi0_mean(bdkp_b0g_mean,dm_addedg2pi0)")
  w.factory("sum::bdkp_bspi0_mean(bdkp_b0g_mean, dm_b02bs, dm_addedg2pi0)")
  # dm shift for johnson and CB comps
  w.factory("dm_j2cb[-140,-200,-80]")
  w.factory("sum::bdkp_b0g_meanj(bdkp_b0g_mean, dm_j2cb)")
  w.factory("sum::bdkp_b0pi0_meanj(bdkp_b0pi0_mean, dm_j2cb)")
  w.factory("sum::bdkp_bsg_meanj(bdkp_bsg_mean, dm_j2cb)")
  w.factory("sum::bdkp_bspi0_meanj(bdkp_bspi0_mean, dm_j2cb)")
  # dm shift for the small part to the big
  w.factory("dm_s2b[-200,-250,-100]")
  w.factory("sum::bdkp_b0g_means(bdkp_b0g_mean, dm_s2b)")
  w.factory("sum::bdkp_b0pi0_means(bdkp_b0pi0_mean, dm_s2b)")
  w.factory("sum::bdkp_bsg_means(bdkp_bsg_mean, dm_s2b)")
  w.factory("sum::bdkp_bspi0_means(bdkp_bspi0_mean, dm_s2b)")

  # sigmas
  w.factory("bdkp_b0g_sigma[30,0,200]")
  w.factory("prod::bdkp_bsg_sigma(bdkp_b0g_sigma, ssig_b02bs)")
  # sig scale for g -> pi0 in this case
  w.factory("ssig_addedg2pi0[0.2,0.,10.]")
  w.factory("prod::bdkp_b0pi0_sigma(bdkp_b0g_sigma, ssig_addedg2pi0)")
  w.factory("prod::bdkp_bspi0_sigma(bdkp_b0g_sigma, ssig_b02bs, ssig_addedg2pi0)")
  # sig scale for small part to the big
  w.factory("ssig_s2b[0.2,0.01,4.]")
  w.factory("prod::bdkp_b0g_sigmas(bdkp_b0g_sigma, ssig_s2b)")
  w.factory("prod::bdkp_b0pi0_sigmas(bdkp_b0pi0_sigma, ssig_s2b)")
  w.factory("prod::bdkp_bsg_sigmas(bdkp_bsg_sigma, ssig_s2b)")
  w.factory("prod::bdkp_bspi0_sigmas(bdkp_bspi0_sigma, ssig_s2b)")

  # tail params

  # Big part is 1CB + Johns :
  w.factory("bdkp_alpha1[2,0,4]")
  w.factory("bdkp_n1[3]")
  w.factory("bdkp_f1[0.4,0.1,1.]")
  w.factory("bdkp_g_gamma[-3,-20,0.]")
  w.factory("bdkp_g_delta[3,1,5]")
  w.factory("bdkp_pi0_gamma[-3,-20,0.]")
  w.factory("bdkp_pi0_delta[3,1,5]")

  # Small part is 1CB
  w.factory("bdkp_alpha2[2,0,4]")
  w.factory("bdkp_n2[3]")
  w.factory("bdkp_f2[0.01,0.,0.1]")

  # B0g signal model
  w.factory("CBShape::bdkp_b0g_cb1( B_DTFDict_D0_B_M, bdkp_b0g_mean, bdkp_b0g_sigma, bdkp_alpha1, bdkp_n1 )")
  w.factory("CBShape::bdkp_b0g_cb2( B_DTFDict_D0_B_M, bdkp_b0g_means, bdkp_b0g_sigmas, bdkp_alpha2, bdkp_n2 )")
  w.factory("Johnson::bdkp_b0g_john( B_DTFDict_D0_B_M, bdkp_b0g_meanj, bdkp_b0g_sigma, bdkp_g_gamma, bdkp_g_delta )")
  #w.factory("SUM::bdkp_mc_pdf_b0g( bdkp_f1*bdkp_b0g_john, bdkp_b0g_cb1 )")
  w.factory("SUM::bdkp_mc_pdf_b0g( bdkp_f1*bdkp_b0g_john, bdkp_f2*bdkp_b0g_cb2, bdkp_b0g_cb1 )")
  # B0pi0 signal model
  w.factory("CBShape::bdkp_b0pi0_cb1( B_DTFDict_D0_B_M, bdkp_b0pi0_mean, bdkp_b0pi0_sigma, bdkp_alpha1, bdkp_n1 )")
  w.factory("CBShape::bdkp_b0pi0_cb2( B_DTFDict_D0_B_M, bdkp_b0pi0_means, bdkp_b0pi0_sigmas, bdkp_alpha2, bdkp_n2 )")
  w.factory("Johnson::bdkp_b0pi0_john( B_DTFDict_D0_B_M, bdkp_b0pi0_meanj, bdkp_b0pi0_sigma, bdkp_pi0_gamma, bdkp_pi0_delta )")
  #w.factory("SUM::bdkp_mc_pdf_b0pi0( bdkp_f1*bdkp_b0pi0_john, bdkp_b0pi0_cb1 )")
  w.factory("SUM::bdkp_mc_pdf_b0pi0( bdkp_f1*bdkp_b0pi0_john, bdkp_f2*bdkp_b0pi0_cb2, bdkp_b0pi0_cb1 )")
  # Bsg signal model
  w.factory("CBShape::bdkp_bsg_cb1( B_DTFDict_D0_B_M, bdkp_bsg_mean, bdkp_bsg_sigma, bdkp_alpha1, bdkp_n1 )")
  w.factory("CBShape::bdkp_bsg_cb2( B_DTFDict_D0_B_M, bdkp_bsg_means, bdkp_bsg_sigmas, bdkp_alpha2, bdkp_n2 )")
  w.factory("Johnson::bdkp_bsg_john( B_DTFDict_D0_B_M, bdkp_bsg_meanj, bdkp_bsg_sigma, bdkp_g_gamma, bdkp_g_delta )")
  #w.factory("SUM::bdkp_mc_pdf_bsg( bdkp_f1*bdkp_bsg_john, bdkp_bsg_cb1 )")
  w.factory("SUM::bdkp_mc_pdf_bsg( bdkp_f1*bdkp_bsg_john, bdkp_f2*bdkp_bsg_cb2, bdkp_bsg_cb1 )")
  # Bspi0 signal model
  w.factory("CBShape::bdkp_bspi0_cb1( B_DTFDict_D0_B_M, bdkp_bspi0_mean, bdkp_bspi0_sigma, bdkp_alpha1, bdkp_n1 )")
  w.factory("CBShape::bdkp_bspi0_cb2( B_DTFDict_D0_B_M, bdkp_bspi0_means, bdkp_bspi0_sigmas, bdkp_alpha2, bdkp_n2 )")
  w.factory("Johnson::bdkp_bspi0_john( B_DTFDict_D0_B_M, bdkp_bspi0_meanj, bdkp_bspi0_sigma, bdkp_pi0_gamma, bdkp_pi0_delta )")
  #w.factory("SUM::bdkp_mc_pdf_bspi0( bdkp_f1*bdkp_bspi0_john, bdkp_bspi0_cb1 )")
  w.factory("SUM::bdkp_mc_pdf_bspi0( bdkp_f1*bdkp_bspi0_john, bdkp_f2*bdkp_bspi0_cb2, bdkp_bspi0_cb1 )")

  # make sim pdf and sim data
  CreateSimPdf(w,'bdkp_mc')
  CreateSimData(w,'bdkp_mc')

  c = TCanvas('c_bdkp_mc','c',1400,1200)
  c.Divide(2,2)

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'bdkp_mc_pdf_%s'%(samp)
      dsetname = 'bdkp_mc_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname))
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp)
    c.Update()
    c.Modified()
    c.Print("plots/bdkp_mc_free_fit.pdf")

  if sim:
    pdfname = 'bdkp_mc_sim_pdf'
    dsetname = 'bdkp_mc_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    for i, samp in enumerate(samples):
      dopts = RooLinkedList()
      dopts.Add( RooCmdArg( rf.Cut('cat==cat::%s'%samp) ) )
      popts = RooLinkedList()
      popts.Add( RooCmdArg( rf.Slice(w.cat("cat"),samp) ) )
      if w.data(dsetname): popts.Add( RooCmdArg( rf.ProjWData(w.data(dsetname)) ) )
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp, popts, dopts)

    c.Update()
    c.Modified()
    c.Print("plots/bdkp_mc_sim_fit.pdf")

    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.saveSnapshot(pdfname.replace('pdf','fit'), pars)

  w.writeToFile('files/w_bdkp.root')

### B -> D*K / D*pi shape ###
### Will be B+ -> D*K+ with a random pi- added for the B0 mode
### or B+ -> D*pi+ with a random K- added for the Bs mode
### So B0 to Bs will have a shift and smear representing a K -> pi addition
### For g -> pi0 the shift and smear should be the same as for the full rec
### signal case
def BDstHShape(free=True,sim=True):
  w = RooWorkspace('w_bdsth','w_bdsth')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))

  files = { 'b0g'   : path+"B0gamma/9_MC_B0gamma_Bg_Bu2DstK_BestCut.root",
            'b0pi0' : path+"B0pi0/9_MC_B0pi0_Bg_Bu2DstK_BestCut.root",
            'bsg'   : path+"Bsgamma/9_MC_Bsgamma_Bg_Bu2Dstpi_BestCut.root",
            'bspi0' : path+"Bspi0/9_MC_Bspi0_Bg_Bu2Dstpi_BestCut.root"
          }

  ### Make the MC dsets ###
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("bdsth_mc_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  ### Make the signal PDF ###
  SharedParams(w)

  # means
  w.factory("bdsth_b0g_mean[5480,5400,5550]")
  w.factory("sum::bdsth_b0pi0_mean(bdsth_b0g_mean,dm_g2pi0)")
  # shift for B0 -> Bs (either missing pi instead of missing K)
  w.factory("dm_addedk2pi[350,250,500]")
  w.factory("sum::bdsth_bsg_mean(bdsth_b0g_mean,dm_addedk2pi)")
  w.factory("sum::bdsth_bspi0_mean(bdsth_b0g_mean,dm_g2pi0,dm_addedk2pi)")

  # sigmas
  w.factory("bdsth_b0g_sigma[30,20,200]")
  w.factory("prod::bdsth_b0pi0_sigma(bdsth_b0g_sigma,ssig_g2pi0)")
  # shift for B0 -> Bs (either missing pi instead of missing K)
  w.factory("ssig_addedk2pi[1.,0.1,5.]")
  w.factory("prod::bdsth_bsg_sigma(bdsth_b0g_sigma,ssig_addedk2pi)")
  w.factory("prod::bdsth_bspi0_sigma(bdsth_b0g_sigma,ssig_g2pi0,ssig_addedk2pi)")

  # tail params

  # Just 1 CB
  w.factory("bdsth_alpha1[-2,-4,0]")
  w.factory("bdsth_n1[3]")
  w.factory("bdsth_alpha2[2,0,4]")
  w.factory("bdsth_n2[3]")
  w.factory("bdsth_f1[0.5,0.,1.]")

  # B0g signal model
  w.factory("CBShape::bdsth_b0g_cb1( B_DTFDict_D0_B_M, bdsth_b0g_mean, bdsth_b0g_sigma, bdsth_alpha1, bdsth_n1 )")
  w.factory("CBShape::bdsth_b0g_cb2( B_DTFDict_D0_B_M, bdsth_b0g_mean, bdsth_b0g_sigma, bdsth_alpha2, bdsth_n2 )")
  w.factory("SUM::bdsth_mc_pdf_b0g( bdsth_f1*bdsth_b0g_cb1, bdsth_b0g_cb2 )")
  # B0pi0 signal model
  w.factory("CBShape::bdsth_b0pi0_cb1( B_DTFDict_D0_B_M, bdsth_b0pi0_mean, bdsth_b0pi0_sigma, bdsth_alpha1, bdsth_n1 )")
  w.factory("CBShape::bdsth_b0pi0_cb2( B_DTFDict_D0_B_M, bdsth_b0pi0_mean, bdsth_b0pi0_sigma, bdsth_alpha2, bdsth_n2 )")
  w.factory("SUM::bdsth_mc_pdf_b0pi0( bdsth_f1*bdsth_b0pi0_cb1, bdsth_b0pi0_cb2 )")
  # Bsg signal model
  w.factory("CBShape::bdsth_bsg_cb1( B_DTFDict_D0_B_M, bdsth_bsg_mean, bdsth_bsg_sigma, bdsth_alpha1, bdsth_n1 )")
  w.factory("CBShape::bdsth_bsg_cb2( B_DTFDict_D0_B_M, bdsth_bsg_mean, bdsth_bsg_sigma, bdsth_alpha2, bdsth_n2 )")
  w.factory("SUM::bdsth_mc_pdf_bsg( bdsth_f1*bdsth_bsg_cb1, bdsth_bsg_cb2 )")
  # Bspi0 signal model
  w.factory("CBShape::bdsth_bspi0_cb1( B_DTFDict_D0_B_M, bdsth_bspi0_mean, bdsth_bspi0_sigma, bdsth_alpha1, bdsth_n1 )")
  w.factory("CBShape::bdsth_bspi0_cb2( B_DTFDict_D0_B_M, bdsth_bspi0_mean, bdsth_bspi0_sigma, bdsth_alpha2, bdsth_n2 )")
  w.factory("SUM::bdsth_mc_pdf_bspi0( bdsth_f1*bdsth_bspi0_cb1, bdsth_bspi0_cb2 )")

  # make sim pdf and sim data
  CreateSimPdf(w,'bdsth_mc')
  CreateSimData(w,'bdsth_mc')

  c = TCanvas('c_bdsth_mc','c',1400,1200)
  c.Divide(2,2)

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'bdsth_mc_pdf_%s'%(samp)
      dsetname = 'bdsth_mc_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname))
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp)
    c.Update()
    c.Modified()
    c.Print("plots/bdsth_mc_free_fit.pdf")

  if sim:
    pdfname = 'bdsth_mc_sim_pdf'
    dsetname = 'bdsth_mc_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    for i, samp in enumerate(samples):
      dopts = RooLinkedList()
      dopts.Add( RooCmdArg( rf.Cut('cat==cat::%s'%samp) ) )
      popts = RooLinkedList()
      popts.Add( RooCmdArg( rf.Slice(w.cat("cat"),samp) ) )
      if w.data(dsetname): popts.Add( RooCmdArg( rf.ProjWData(w.data(dsetname)) ) )
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp, popts, dopts)

    c.Update()
    c.Modified()
    c.Print("plots/bdsth_mc_sim_fit.pdf")

    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.saveSnapshot(pdfname.replace('pdf','fit'), pars)

  w.writeToFile('files/w_bdsth.root')

### Part-Rec Shape ###

def PartRecShape(free=True,sim=True):
  w = RooWorkspace('w_partrec','w_partrec')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))

  files = { 'b0g'   : path+"B0gamma/Bu2DstK1_tree_MeV_FitRegion_Small.root",
            'b0pi0' : path+"B0pi0/Bu2DstK1_tree_MeV_FitRegion_Small.root",
            'bsg'   : path+"Bsgamma/Bs2DstK1_tree_MeV_InRange_Small.root",
            'bspi0' : path+"Bspi0/Bs2DstK1_tree_MeV_InRange_Small.root"
          }

  ### Make the MC dsets ###
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("partrec_mc_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  ### Make the signal PDF ###
  SharedParams(w)

  # means
  w.factory("partrec_b0_m0[5260,5150,5300]")
  w.factory("sum::partrec_bs_m0( partrec_b0_m0, dm_b02bs )")
  w.factory("partrec_c[10,0.1,100]")
  w.factory("partrec_p[1,0,100]")

  # B0g signal model
  #w.factory("ArgusBG::partrec_mc_pdf_b0g( partrec_m, partrec_m0, partrec_p, partrec_c  )")
  # B0pi0 signal model
  #w.factory("ArgusBG::partrec_mc_pdf_b0pi0( partrec_m, partrec_m0, partrec_p, partrec_c )")
  # Bsg signal model
  #w.factory("ArgusBG::partrec_mc_pdf_bsg( partrec_m, partrec_m0, partrec_p, partrec_c )")
  # Bspi0 signal model
  #w.factory("ArgusBG::partrec_mc_pdf_bspi0( partrec_m, partrec_m0, partrec_p, partrec_c )")

  # means
  w.factory("partrec_b0g_mean[4900,4500,5100]")
  w.factory("sum::partrec_b0pi0_mean(partrec_b0g_mean,dm_g2pi0)")
  w.factory("sum::partrec_bsg_mean(partrec_b0g_mean,dm_b02bs)")
  w.factory("sum::partrec_bspi0_mean(partrec_b0g_mean,dm_b02bs,dm_g2pi0)")

  # sigmas
  w.factory("partrec_b0g_sigma[30,20,200]")
  w.factory("prod::partrec_b0pi0_sigma(partrec_b0g_sigma,ssig_g2pi0)")
  w.factory("prod::partrec_bsg_sigma(partrec_b0g_sigma,ssig_b02bs)")
  w.factory("prod::partrec_bspi0_sigma(partrec_b0g_sigma,ssig_g2pi0,ssig_b02bs)")
  w.factory("partrec_b0g_sigma2[30,20,200]")
  w.factory("prod::partrec_b0pi0_sigma2(partrec_b0g_sigma2,ssig_g2pi0)")
  w.factory("prod::partrec_bsg_sigma2(partrec_b0g_sigma2,ssig_b02bs)")
  w.factory("prod::partrec_bspi0_sigma2(partrec_b0g_sigma2,ssig_g2pi0,ssig_b02bs)")

  # tail params

  # Just 1 CB
  w.factory("partrec_alpha1[-2,-4,0]")
  w.factory("partrec_n1[3]")
  w.factory("partrec_alpha2[2,0,4]")
  w.factory("partrec_n2[3]")
  w.factory("partrec_f1[0.5,0.,1.]")

  # B0g signal model
  #w.factory("CBShape::partrec_b0g_cb1( B_DTFDict_D0_B_M, partrec_b0g_mean, partrec_b0g_sigma, partrec_alpha1, partrec_n1 )")
  w.factory("CBShape::partrec_b0g_cb2( B_DTFDict_D0_B_M, partrec_b0g_mean, partrec_b0g_sigma, partrec_alpha2, partrec_n2 )")
  #w.factory("Gaussian::partrec_b0g_gaus( B_DTFDict_D0_B_M, partrec_b0g_mean, partrec_b0g_sigma2 )")
  w.factory("ArgusBG::partrec_b0g_bg( B_DTFDict_D0_B_M, partrec_b0_m0, partrec_c, partrec_p )")
  w.factory("SUM::partrec_mc_pdf_b0g( partrec_f1*partrec_b0g_cb2, partrec_b0g_bg )")
  # B0pi0 signal model
  #w.factory("CBShape::partrec_b0pi0_cb1( B_DTFDict_D0_B_M, partrec_b0pi0_mean, partrec_b0pi0_sigma, partrec_alpha1, partrec_n1 )")
  w.factory("CBShape::partrec_b0pi0_cb2( B_DTFDict_D0_B_M, partrec_b0pi0_mean, partrec_b0pi0_sigma, partrec_alpha2, partrec_n2 )")
  #w.factory("Gaussian::partrec_b0pi0_gaus( B_DTFDict_D0_B_M, partrec_b0pi0_mean, partrec_b0pi0_sigma2 )")
  w.factory("ArgusBG::partrec_b0pi0_bg( B_DTFDict_D0_B_M, partrec_b0_m0, partrec_c, partrec_p )")
  w.factory("SUM::partrec_mc_pdf_b0pi0( partrec_f1*partrec_b0pi0_cb2, partrec_b0pi0_bg )")
  # Bsg signal model
  #w.factory("CBShape::partrec_bsg_cb1( B_DTFDict_D0_B_M, partrec_bsg_mean, partrec_bsg_sigma, partrec_alpha1, partrec_n1 )")
  w.factory("CBShape::partrec_bsg_cb2( B_DTFDict_D0_B_M, partrec_bsg_mean, partrec_bsg_sigma, partrec_alpha2, partrec_n2 )")
  #w.factory("Gaussian::partrec_bsg_gaus( B_DTFDict_D0_B_M, partrec_bsg_mean, partrec_bsg_sigma2 )")
  w.factory("ArgusBG::partrec_bsg_bg( B_DTFDict_D0_B_M, partrec_bs_m0, partrec_c, partrec_p )")
  w.factory("SUM::partrec_mc_pdf_bsg( partrec_f1*partrec_bsg_cb2, partrec_bsg_bg )")
  # Bspi0 signal model
  #w.factory("CBShape::partrec_bspi0_cb1( B_DTFDict_D0_B_M, partrec_bspi0_mean, partrec_bspi0_sigma, partrec_alpha1, partrec_n1 )")
  w.factory("CBShape::partrec_bspi0_cb2( B_DTFDict_D0_B_M, partrec_bspi0_mean, partrec_bspi0_sigma, partrec_alpha2, partrec_n2 )")
  #w.factory("Gaussian::partrec_bspi0_gaus( B_DTFDict_D0_B_M, partrec_bspi0_mean, partrec_bspi0_sigma2 )")
  w.factory("ArgusBG::partrec_bspi0_bg( B_DTFDict_D0_B_M, partrec_bs_m0, partrec_c, partrec_p )")
  w.factory("SUM::partrec_mc_pdf_bspi0( partrec_f1*partrec_bspi0_cb2, partrec_bspi0_bg )")

  # make sim pdf and sim data
  CreateSimPdf(w,'partrec_mc')
  CreateSimData(w,'partrec_mc')

  c = TCanvas('c_partrec_mc','c',1400,1200)
  c.Divide(2,2)

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'partrec_mc_pdf_%s'%(samp)
      dsetname = 'partrec_mc_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname))
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp)
    c.Update()
    c.Modified()
    c.Print("plots/partrec_mc_free_fit.pdf")

  if sim:
    pdfname = 'partrec_mc_sim_pdf'
    dsetname = 'partrec_mc_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    for i, samp in enumerate(samples):
      dopts = RooLinkedList()
      dopts.Add( RooCmdArg( rf.Cut('cat==cat::%s'%samp) ) )
      popts = RooLinkedList()
      popts.Add( RooCmdArg( rf.Slice(w.cat("cat"),samp) ) )
      if w.data(dsetname): popts.Add( RooCmdArg( rf.ProjWData(w.data(dsetname)) ) )
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp, popts, dopts)

    c.Update()
    c.Modified()
    c.Print("plots/partrec_mc_sim_fit.pdf")

    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.saveSnapshot(pdfname.replace('pdf','fit'), pars)

  w.writeToFile('files/w_partrec.root')

### Lb -> D0ph ###
### h will be pi- for B0 and K- for Bs
### so for B0 we have proton misID as kaon and Bs we have proton misID as pi
### and then add a random g or pi0
### suspect this maybe mixed up with the D* sample?
def Lb2DPHShape(free=True,sim=True):
  w = RooWorkspace('w_lbdph','w_lbdph')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))

  files = { 'b0g'   : path+"B0gamma/9_MC_B0gamma_Bg_Lb2DpK_BestCut.root",
            'b0pi0' : path+"B0pi0/9_MC_B0pi0_Bg_Lb2DpK_BestCut.root",
            'bsg'   : path+"Bsgamma/9_MC_Bsgamma_Bg_Lb2Dppi_BestCut.root",
            'bspi0' : path+"Bspi0/9_MC_Bspi0_Bg_Lb2Dppi_BestCut.root"
          }

  ### Make the MC dsets ###
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("lbdph_mc_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  ### Make the signal PDF ###
  SharedParams(w)

  # means
  w.factory("lbdph_b0g_mean[5480,5400,5550]")
  w.factory("sum::lbdph_b0pi0_mean(lbdph_b0g_mean,dm_g2pi0)")
  # shift for B0 -> Bs (either misID pi for proton of misID K for proton)
  w.factory("dm_mispk2pi[100,20,200]")
  w.factory("sum::lbdph_bsg_mean(lbdph_b0g_mean,dm_mispk2pi)")
  w.factory("sum::lbdph_bspi0_mean(lbdph_b0g_mean,dm_g2pi0,dm_mispk2pi)")

  # sigmas
  w.factory("lbdph_b0g_sigma[30,20,500]")
  w.factory("prod::lbdph_b0pi0_sigma(lbdph_b0g_sigma,ssig_g2pi0)")
  # shift for B0 -> Bs (either missing pi instead of missing K)
  w.factory("ssig_mispk2p[1.,0.1,5.]")
  w.factory("prod::lbdph_bsg_sigma(lbdph_b0g_sigma,ssig_mispk2p)")
  w.factory("prod::lbdph_bspi0_sigma(lbdph_b0g_sigma,ssig_g2pi0,ssig_mispk2p)")

  # tail params

  # Just 1 CB
  w.factory("lbdph_alpha1[-2,-4,0]")
  w.factory("lbdph_n1[3]")
  w.factory("lbdph_alpha2[2,0,4]")
  w.factory("lbdph_n2[3]")
  w.factory("lbdph_f1[0.5,0.,1.]")

  # B0g signal model
  w.factory("CBShape::lbdph_b0g_cb1( B_DTFDict_D0_B_M, lbdph_b0g_mean, lbdph_b0g_sigma, lbdph_alpha1, lbdph_n1 )")
  w.factory("CBShape::lbdph_b0g_cb2( B_DTFDict_D0_B_M, lbdph_b0g_mean, lbdph_b0g_sigma, lbdph_alpha2, lbdph_n2 )")
  w.factory("SUM::lbdph_mc_pdf_b0g( lbdph_f1*lbdph_b0g_cb1, lbdph_b0g_cb2 )")
  # B0pi0 signal model
  w.factory("CBShape::lbdph_b0pi0_cb1( B_DTFDict_D0_B_M, lbdph_b0pi0_mean, lbdph_b0pi0_sigma, lbdph_alpha1, lbdph_n1 )")
  w.factory("CBShape::lbdph_b0pi0_cb2( B_DTFDict_D0_B_M, lbdph_b0pi0_mean, lbdph_b0pi0_sigma, lbdph_alpha2, lbdph_n2 )")
  w.factory("SUM::lbdph_mc_pdf_b0pi0( lbdph_f1*lbdph_b0pi0_cb1, lbdph_b0pi0_cb2 )")
  # Bsg signal model
  w.factory("CBShape::lbdph_bsg_cb1( B_DTFDict_D0_B_M, lbdph_bsg_mean, lbdph_bsg_sigma, lbdph_alpha1, lbdph_n1 )")
  w.factory("CBShape::lbdph_bsg_cb2( B_DTFDict_D0_B_M, lbdph_bsg_mean, lbdph_bsg_sigma, lbdph_alpha2, lbdph_n2 )")
  w.factory("SUM::lbdph_mc_pdf_bsg( lbdph_f1*lbdph_bsg_cb1, lbdph_bsg_cb2 )")
  # Bspi0 signal model
  w.factory("CBShape::lbdph_bspi0_cb1( B_DTFDict_D0_B_M, lbdph_bspi0_mean, lbdph_bspi0_sigma, lbdph_alpha1, lbdph_n1 )")
  w.factory("CBShape::lbdph_bspi0_cb2( B_DTFDict_D0_B_M, lbdph_bspi0_mean, lbdph_bspi0_sigma, lbdph_alpha2, lbdph_n2 )")
  w.factory("SUM::lbdph_mc_pdf_bspi0( lbdph_f1*lbdph_bspi0_cb1, lbdph_bspi0_cb2 )")

  # make sim pdf and sim data
  CreateSimPdf(w,'lbdph_mc')
  CreateSimData(w,'lbdph_mc')

  c = TCanvas('c_lbdph_mc','c',1400,1200)
  c.Divide(2,2)

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'lbdph_mc_pdf_%s'%(samp)
      dsetname = 'lbdph_mc_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname))
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp)
    c.Update()
    c.Modified()
    c.Print("plots/lbdph_mc_free_fit.pdf")

  if sim:
    pdfname = 'lbdph_mc_sim_pdf'
    dsetname = 'lbdph_mc_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    for i, samp in enumerate(samples):
      dopts = RooLinkedList()
      dopts.Add( RooCmdArg( rf.Cut('cat==cat::%s'%samp) ) )
      popts = RooLinkedList()
      popts.Add( RooCmdArg( rf.Slice(w.cat("cat"),samp) ) )
      if w.data(dsetname): popts.Add( RooCmdArg( rf.ProjWData(w.data(dsetname)) ) )
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp, popts, dopts)

    c.Update()
    c.Modified()
    c.Print("plots/lbdph_mc_sim_fit.pdf")

    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.saveSnapshot(pdfname.replace('pdf','fit'), pars)

  w.writeToFile('files/w_lbdph.root')

### Lb -> D*ph ###
### h will be pi- for B0 and K- for Bs
### so for B0 we have proton misID as kaon and Bs we have proton misID as pi
def Lb2DstPHShape(free=True,sim=True):
  w = RooWorkspace('w_lbdstph','w_lbdstph')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))

  files = { 'b0g'   : path+"B0gamma/9_MC_B0gamma_Bg_Lb2DstpK_BestCut.root",
            'b0pi0' : path+"B0pi0/9_MC_B0pi0_Bg_Lb2DstpK_BestCut.root",
            'bsg'   : path+"Bsgamma/9_MC_Bsgamma_Bg_Lb2Dstppi_BestCut.root",
            'bspi0' : path+"Bspi0/9_MC_Bspi0_Bg_Lb2Dstppi_BestCut.root"
          }

  ### Make the MC dsets ###
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("lbdstph_mc_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  ### Make the signal PDF ###
  # going to make life easy and RooKeys them
  for i, samp in enumerate(samples):
    pdfname  = 'lbdstph_mc_pdf_%s'%(samp)
    dsetname = 'lbdstph_mc_%s'%(samp)
    pdf = RooKeysPdf(pdfname,pdfname,w.var("B_DTFDict_D0_B_M"),w.data(dsetname))
    getattr(w,'import')(pdf)


  # make sim pdf and sim data
  CreateSimPdf(w,'lbdstph_mc')
  CreateSimData(w,'lbdstph_mc')

  c = TCanvas('c_lbdstph_mc','c',1400,1200)
  c.Divide(2,2)

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'lbdstph_mc_pdf_%s'%(samp)
      dsetname = 'lbdstph_mc_%s'%(samp)
      #w.pdf(pdfname).fitTo(w.data(dsetname)) # nothing to fit in this case
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp)
    c.Update()
    c.Modified()
    c.Print("plots/lbdstph_mc_free_fit.pdf")

  if sim:
    pdfname = 'lbdstph_mc_sim_pdf'
    dsetname = 'lbdstph_mc_sim_data'
    #w.pdf(pdfname).fitTo( w.data(dsetname) ) # nothing to fit in this case
    for i, samp in enumerate(samples):
      dopts = RooLinkedList()
      dopts.Add( RooCmdArg( rf.Cut('cat==cat::%s'%samp) ) )
      popts = RooLinkedList()
      popts.Add( RooCmdArg( rf.Slice(w.cat("cat"),samp) ) )
      if w.data(dsetname): popts.Add( RooCmdArg( rf.ProjWData(w.data(dsetname)) ) )
      plot( c.cd(i+1), w, pdfname, dsetname, True, samp, popts, dopts)

    c.Update()
    c.Modified()
    c.Print("plots/lbdstph_mc_sim_fit.pdf")

    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.saveSnapshot(pdfname.replace('pdf','fit'), pars)

  w.writeToFile('files/w_lbdstph.root')

