import os
from ROOT import RooFit as rf
from ROOT import TFile, TTree, TCanvas
from ROOT import RooWorkspace, RooCmdArg
from utils import *

path = "/Users/matt/Scratch/lhcb/ArnauFits/Tuples/"

def ImportMCShapes(w):

  mc_types = [ 'sig', 'misrec', 'bdstpp', 'bdstkk', 'bdkp', 'bdsth', 'partrec', 'lbdph', 'lbdstph' ]

  for typ in mc_types:
    assert( os.path.exists( 'files/w_%s.root'%typ ) )
    tf = TFile( 'files/w_%s.root'%typ )
    wmc = tf.Get('w_%s'%typ)
    assert(wmc)
    wmc.loadSnapshot('%s_mc_sim_fit'%typ)
    pdfname = '%s_mc_sim_pdf'%typ
    getattr(w,'import')(wmc.pdf(pdfname), rf.RecycleConflictNodes() )
    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.defineSet(pdfname+'_pars',pars)
    tf.Close()

  # some special jiggery pokery for the dstkk stuff (as it was fitted in 8 categories not 4 to MC)
  for samp in [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]:
    bdstkk = w.pdf('bdstkk_mc_pdf_b0_%s'%samp)
    bdstkk.SetName('bdstkk_mc_pdf_%s'%samp)
    bsdstkk = w.pdf('bdstkk_mc_pdf_bs_%s'%samp)
    bsdstkk.SetName('bsdstkk_mc_pdf_%s'%samp)
    getattr(w,'import')(bdstkk, rf.RecycleConflictNodes() )
    bd_pars = w.pdf('bdstkk_mc_pdf_%s'%samp).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.defineSet('bdstkk_mc_pdf_%s_pars'%samp,bd_pars)
    getattr(w,'import')(bsdstkk, rf.RecycleConflictNodes() )
    bs_pars = w.pdf('bsdstkk_mc_pdf_%s'%samp).getParameters(RooArgSet(w.var('B_DTFDict_D0_B_M')))
    w.defineSet('bsdstkk_mc_pdf_%s_pars'%samp,bs_pars)

def DataFit(free=True, sim=True):

  w = RooWorkspace('w_data','w_data')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))


  files = { 'b0g'   : path+"B0gamma/9_B2DstKpi_Dst2DgammTuple_BestCut.root",
            'b0pi0' : path+"B0pi0/9_B2Dstpi0Tuple_BestCut.root",
            'bsg'   : path+"Bsgamma/9_Bs2DstKpi_Dst2DgammaTuple_BestCut.root",
            'bspi0' : path+"Bspi0/9_Bs2Dstpi0Tuple_BestCut.root"
          }

  # Make the dsets
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("data_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  # Make the total pdf

  # First try and merge the different bits from the different workspaces
  ImportMCShapes(w)

  # Then make combinatorial shape in each cateogry
  # let these be independent for now
  for samp in samples:
    w.factory("comb_mc_%s_p0[-0.001,-0.1,0.]"%samp)
    w.factory("Exponential::comb_mc_pdf_%s( B_DTFDict_D0_B_M, comb_mc_%s_p0 )"%(samp,samp) )
    w.factory("%s_comb_y[3000,0,12000]"%samp)

  # Now need to figure out what yields to restrict

  # sig yield first (require b0 / bs ratio consistent between g and pi0)
  w.factory("b0g_sig_y[3000,0,12000]")
  w.factory("b0pi0_sig_y[800,0,4000]")
  w.factory("bs2b0_rat[2.5,1.,4.]")
  w.factory("prod::bsg_sig_y(b0g_sig_y, bs2b0_rat)")
  w.factory("prod::bspi0_sig_y(b0pi0_sig_y, bs2b0_rat)")

  # now mis rec yield (ratio of this to sig should be the same for b0 and bs but will be different for g vs pi0)
  w.factory("misrec_to_sig_rat_g[0.2,0.001,0.6]")
  w.factory("misrec_to_sig_rat_pi0[0.2,0.001,0.6]")
  w.factory("prod::b0g_misrec_y( misrec_to_sig_rat_g, b0g_sig_y )")
  w.factory("prod::b0pi0_misrec_y( misrec_to_sig_rat_pi0, b0pi0_sig_y )")
  w.factory("prod::bsg_misrec_y( misrec_to_sig_rat_g, bsg_sig_y )")
  w.factory("prod::bspi0_misrec_y( misrec_to_sig_rat_pi0, bspi0_sig_y )")

  # the cases of B->D*pipi, B->D*KK, Lb->D*ph all involve a misID so will
  # be different for B0 and Bs (as they differ with a K or pi misID) however
  # for all of these the ratio of g -> pi0 should be the same
  # there is also Bs->D*KK which should scale the same for g and pi0 modes
  w.factory("misid_g2pi0_rat[0.1,0.0001,10.]")
  w.factory("b0g_bdstpp_y[1000,0,12000]")
  w.factory("bsg_bdstpp_y[1000,0,12000]")
  w.factory("prod::b0pi0_bdstpp_y( misid_g2pi0_rat, b0g_bdstpp_y )")
  w.factory("prod::bspi0_bdstpp_y( misid_g2pi0_rat, bsg_bdstpp_y )")
  w.factory("b0g_bdstkk_y[1000,0,12000]")
  w.factory("bsg_bdstkk_y[1000,0,12000]")
  w.factory("prod::b0pi0_bdstkk_y( misid_g2pi0_rat, b0g_bdstkk_y )")
  w.factory("prod::bspi0_bdstkk_y( misid_g2pi0_rat, bsg_bdstkk_y )")
  w.factory("b0g_lbdstph_y[1000,0,12000]")
  w.factory("bsg_lbdstph_y[1000,0,12000]")
  w.factory("prod::b0pi0_lbdstph_y( misid_g2pi0_rat, b0g_lbdstph_y )")
  w.factory("prod::bspi0_lbdstph_y( misid_g2pi0_rat, bsg_lbdstph_y )")
  w.factory("bsdstkk_to_bdstkk_rat[1.,0.1,2.]")
  w.factory("prod::b0g_bsdstkk_y( bsdstkk_to_bdstkk_rat, b0g_bdstkk_y )")
  w.factory("prod::b0pi0_bsdstkk_y( bsdstkk_to_bdstkk_rat, b0pi0_bdstkk_y )")
  w.factory("prod::bsg_bsdstkk_y( bsdstkk_to_bdstkk_rat, bsg_bdstkk_y )")
  w.factory("prod::bspi0_bsdstkk_y( bsdstkk_to_bdstkk_rat, bspi0_bdstkk_y )")

  # B -> DKpi same logic as misrec
  w.factory("bdkp_to_sig_rat_g[0.2,0.001,0.6]")
  w.factory("bdkp_to_sig_rat_pi0[0.2,0.001,0.6]")
  w.factory("prod::b0g_bdkp_y( bdkp_to_sig_rat_g, b0g_sig_y )")
  w.factory("prod::b0pi0_bdkp_y( bdkp_to_sig_rat_pi0, b0pi0_sig_y )")
  w.factory("prod::bsg_bdkp_y( bdkp_to_sig_rat_g, bsg_sig_y )")
  w.factory("prod::bspi0_bdkp_y( bdkp_to_sig_rat_pi0, bspi0_sig_y )")

  # B -> D* K / B -> D* pi (adding random pi- to B0 and random K- to Bs0)
  # so ratio to signal should be same for both g and pi0 modes but
  # different for B0 -> Bs
  w.factory("bdsth_to_sig_rat_addpi[0.2,0.001,0.6]")
  w.factory("bdsth_to_sig_rat_addk[0.2,0.001,0.6]")
  w.factory("prod::b0g_bdsth_y( bdsth_to_sig_rat_addpi, b0g_sig_y )")
  w.factory("prod::b0pi0_bdsth_y( bdsth_to_sig_rat_addpi, b0pi0_sig_y )")
  w.factory("prod::bsg_bdsth_y( bdsth_to_sig_rat_addk, bsg_sig_y )")
  w.factory("prod::bspi0_bdsth_y( bdsth_to_sig_rat_addk, bspi0_sig_y )")

  # Lb -> Dph (mid-ID k for p and pi for p and add random g or pi0) so will be different for all 4 really
  # express this ratio to the Lb -> D*ph one (they should be similar in magnitude?)
  w.factory("lbdph_to_lbdstph_b0g[1.,0.5,2.]")
  w.factory("lbdph_to_lbdstph_b0pi0[1.,0.5,2.]")
  w.factory("lbdph_to_lbdstph_bsg[1.,0.5,2.]")
  w.factory("lbdph_to_lbdstph_bspi0[1.,0.5,2.]")
  w.factory("prod::b0g_lbdph_y( lbdph_to_lbdstph_b0g, b0g_lbdstph_y )")
  w.factory("prod::b0pi0_lbdph_y( lbdph_to_lbdstph_b0pi0, b0pi0_lbdstph_y )")
  w.factory("prod::bsg_lbdph_y( lbdph_to_lbdstph_bsg, bsg_lbdstph_y )")
  w.factory("prod::bspi0_lbdph_y( lbdph_to_lbdstph_bspi0, bspi0_lbdstph_y )")

  # Part reco shape should have same Bs / B0 ratio
  w.factory("partrec_to_sig_rat_g[0.2,0.001,0.6]")
  w.factory("partrec_to_sig_rat_pi0[0.2,0.001,0.6]")
  w.factory("prod::b0g_partrec_y( partrec_to_sig_rat_g, b0g_sig_y )")
  w.factory("prod::b0pi0_partrec_y( partrec_to_sig_rat_pi0, b0pi0_sig_y )")
  w.factory("prod::bsg_partrec_y( partrec_to_sig_rat_g, bsg_sig_y )")
  w.factory("prod::bspi0_partrec_y( partrec_to_sig_rat_pi0, bspi0_sig_y )")

  components = [ 'sig', 'misrec', 'bdstpp', 'bdstkk', 'bsdstkk', 'bdkp', 'bdsth', 'partrec', 'lbdph', 'lbdstph', 'comb' ]

  for samp in samples:
    fact_str = "SUM::data_pdf_%s("%samp
    for comp in components:
      fact_str += "%s_%s_y*%s_mc_pdf_%s,"%(samp,comp,comp,samp)
    fact_str = fact_str[:-1] +")"
    w.factory(fact_str)
    w.pdf('data_pdf_%s'%samp).Print('v')

  CreateSimPdf(w,'data')
  CreateSimData(w,'data')

  # Now fix appropriate parameters
  # To start with we'll fix all shape parameters from MC and just float the yields (and exponential slope)
  for comp in components:
    if comp == 'comb': continue # no pre-defined shape for combinatorial
    if comp == 'bsdstkk': continue # this params for this piece are covered by bdstkk
    w.set('%s_mc_sim_pdf_pars'%comp).setAttribAll("Constant")

  # Now relax the constraints on a few important params
  w.var("b0g_mean").setConstant(False)
  w.var("b0g_sigma").setConstant(False)
  #w.var("dm_b02bs").setConstant(False)
  #w.var("dm_g2pi0").setConstant(False)
  #w.var("ssig_b02bs").setConstant(False)
  #w.var("ssig_g2pi0").setConstant(False)

  #w.var("b0g_misrec_mean").setConstant(False)
  #w.var("b0g_misrec_sigma").setConstant(False)
  #w.var("dm_missg2addg").setConstant(False)
  #w.var("ssig_missg2addg").setConstant(False)

  w.Print('v')
  w.pdf('data_sim_pdf').Print('v')
  w.data('data_sim_data').Print('v')

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'data_pdf_%s'%(samp)
      dsetname = 'data_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname)) # nothing to fit in this case
      pars = w.pdf(pdfname).getParameters(RooArgSet(w.var("B_DTFDict_D0_B_M")))
      w.saveSnapshot('data_free_fit_%s'%samp,pars)

  if sim:
    pdfname = 'data_sim_pdf'
    dsetname = 'data_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var("B_DTFDict_D0_B_M")))
    w.saveSnapshot('data_sim_fit',pars)

  w.writeToFile('files/w_data.root')

def DataPlot(free=True, sim=True):

  tf = TFile('files/w_data.root')
  w = tf.Get('w_data')

  samples = [ 'b0g', 'b0pi0', 'bsg', 'bspi0' ]

  cs = []

  # free fit
  if free:
    for i, samp in enumerate(samples):
      cs.append( TCanvas('c_data_free_%s'%samp,'',700,600) )
      pdfname  = 'data_pdf_%s'%(samp)
      dsetname = 'data_%s'%(samp)
      w.loadSnapshot('data_free_fit_%s'%samp)
      plotData( cs[-1], w, pdfname, dsetname, True, samp)
      cs[-1].Update()
      cs[-1].Modified()
      cs[-1].Print("plots/data_free_fit_%s.pdf"%samp)
      plotData( cs[-1], w, pdfname, dsetname, True, samp, True)
      cs[-1].Update()
      cs[-1].Modified()
      cs[-1].Print("plots/data_free_fit_%s_log.pdf"%samp)

  # sim fit
  if sim:
    w.loadSnapshot('data_sim_fit')
    for i, samp in enumerate(samples):
      cs.append( TCanvas('c_data_sim_%s'%samp,'',700,600) )
      pdfname  = 'data_pdf_%s'%(samp)
      dsetname = 'data_%s'%(samp)
      plotData( cs[-1], w, pdfname, dsetname, True, samp)
      cs[-1].Update()
      cs[-1].Modified()
      cs[-1].Print("plots/data_sim_fit_%s.pdf"%samp)
      plotData( cs[-1], w, pdfname, dsetname, True, samp, True)
      cs[-1].Update()
      cs[-1].Modified()
      cs[-1].Print("plots/data_sim_fit_%s_log.pdf"%samp)

def ControlDataFit(free=True, sim=True):

  w = RooWorkspace('w_ctrl','w_ctrl')

  samples = [ 'b0g', 'bsg', 'cg', 'b0pi0', 'bspi0', 'cpi0' ]

  # create the category
  w.factory('cat[%s]'%(','.join(samples)))


  files = { 'b0g'   : path+"../New/Tuples/Data/9_B2DstKpi_Dst2DgammTuple_BestCut.root",
            'b0pi0' : path+"../New/Tuples/Data/9_B2Dstpi0Tuple_BestCut.root",
            'bsg'   : path+"../New/Tuples/Data/9_Bs2DstKpi_Dst2DgammaTuple_BestCut.root",
            'bspi0' : path+"../New/Tuples/Data/9_Bs2Dstpi0Tuple_BestCut.root",
            'cg'    : path+"../New/Tuples/Data/9_B2Dstpipi_Dst2DgammTuple_BestCut.root",
            'cpi0'  : path+"../New/Tuples/Data/9_B2Dstpipi_Dst2Dpi0Tuple_No16_BestCut.root"
          }

  # Make the dsets
  w.factory("B_DTFDict_D0_B_M[5100,5900]")
  w.var("B_DTFDict_D0_B_M").setBins(80)
  for samp in samples:
    assert( os.path.exists(files[samp]) )
    tf = TFile( files[samp] )
    t = tf.Get('DecayTree')
    t.SetBranchStatus("*",0)
    t.SetBranchStatus("B_DTFDict_D0_B_M",1)
    dset = RooDataSet("data_%s"%(samp),"",t,RooArgSet(w.var("B_DTFDict_D0_B_M")))
    getattr(w,'import')(dset)
    tf.Close()

  c = TCanvas('ctrl','ctrl',2100,1200)
  c.Divide(3,2)
  for i, samp in enumerate(samples):
    pdfname  = 'data_pdf_%s'%(samp)
    dsetname = 'data_%s'%(samp)
    plot( c.cd(i+1), w, pdfname, dsetname )
  c.Update()
  c.Modified()
  c.Print("plots/ctrl.pdf")

  # Make the total pdf

  # First try and merge the different bits from the different workspaces
  ImportMCShapes(w)

  # now want to make the ctrl pdfs
  # signal and misrec are the same as b0
  sig_cg   = w.pdf('sig_mc_pdf_b0g').Clone('sig_mc_pdf_cg')
  getattr(w,'import')(sig_cg)
  sig_cpi0 = w.pdf('sig_mc_pdf_b0pi0').Clone('sig_mc_pdf_cpi0')
  getattr(w,'import')(sig_cpi0)
  misrec_cg   = w.pdf('misrec_mc_pdf_b0g').Clone('misrec_mc_pdf_cg')
  getattr(w,'import')(misrec_cg)
  misrec_cpi0 = w.pdf('misrec_mc_pdf_b0pi0').Clone('misrec_mc_pdf_cpi0')
  getattr(w,'import')(misrec_cpi0)
  # ignore the lambdas for now
  # there will be some Bs0 -> D*0 piK stuff with a misID'd K (I think this is the big bit which appears below the peak)
  # do 1CB for this
  w.factory('dm_ctrl_misidk2pi[-150,-450,-50]')
  w.factory('sum::bdstkp_cg_mean( b0g_mean, dm_ctrl_misidk2pi)')
  w.factory('sum::bdstkp_cpi0_mean( b0pi0_mean, dm_ctrl_misidk2pi)')
  w.factory('bdstkp_cg_sigma[40,5,200]')
  w.factory('bdstkp_cpi0_sigma[40,5,200]')
  w.factory('bdstkp_c_alpha1[-2.1,-4,0]')
  w.factory('bdstkp_c_n1[3]')
  w.factory('CBShape::bdstkp_mc_pdf_cg( B_DTFDict_D0_B_M, bdstkp_cg_mean, bdstkp_cg_sigma, bdstkp_c_alpha1, bdstkp_c_n1 )')
  w.factory('CBShape::bdstkp_mc_pdf_cpi0( B_DTFDict_D0_B_M, bdstkp_cpi0_mean, bdstkp_cpi0_sigma, bdstkp_c_alpha1, bdstkp_c_n1 )')

  # the stuff above the peak is probably mostly B- -> D*0 pi- with another random pi so could get this shape and shift it
  w.factory('dm_ctrl_dstpi[100,10,300]')
  w.factory('sum::bdsth_cg_mean( bdsth_b0g_mean, dm_ctrl_dstpi )')
  w.factory('sum::bdsth_cpi0_mean( bdsth_b0pi0_mean, dm_ctrl_dstpi )')
  w.factory("CBShape::bdsth_cg_cb1( B_DTFDict_D0_B_M, bdsth_cg_mean, bdsth_b0g_sigma, bdsth_alpha1, bdsth_n1 )")
  w.factory("CBShape::bdsth_cg_cb2( B_DTFDict_D0_B_M, bdsth_cg_mean, bdsth_b0g_sigma, bdsth_alpha2, bdsth_n2 )")
  w.factory("SUM::bdsth_mc_pdf_cg( bdsth_f1*bdsth_cg_cb1, bdsth_cg_cb2 )")
  w.factory("CBShape::bdsth_cpi0_cb1( B_DTFDict_D0_B_M, bdsth_cpi0_mean, bdsth_b0pi0_sigma, bdsth_alpha1, bdsth_n1 )")
  w.factory("CBShape::bdsth_cpi0_cb2( B_DTFDict_D0_B_M, bdsth_cpi0_mean, bdsth_b0pi0_sigma, bdsth_alpha2, bdsth_n2 )")
  w.factory("SUM::bdsth_mc_pdf_cpi0( bdsth_f1*bdsth_cpi0_cb1, bdsth_cpi0_cb2 )")

  # Then make combinatorial shape in each cateogry
  # let these be independent for now
  for samp in samples:
    w.factory("comb_mc_%s_p0[-0.001,-0.1,0.]"%samp)
    w.factory("Exponential::comb_mc_pdf_%s( B_DTFDict_D0_B_M, comb_mc_%s_p0 )"%(samp,samp) )
    w.factory("%s_comb_y[3000,0,12000]"%samp)

  # Now need to figure out what yields to restrict

  # sig yield first (require b0 / bs ratio consistent between g and pi0)
  w.factory("b0g_sig_y[3000,0,12000]")
  w.factory("b0pi0_sig_y[800,0,4000]")
  w.factory("bs2b0_rat[2.5,1.,4.]")
  w.factory("prod::bsg_sig_y(b0g_sig_y, bs2b0_rat)")
  w.factory("prod::bspi0_sig_y(b0pi0_sig_y, bs2b0_rat)")
  w.factory("cg_sig_y[6000,0,32000]")
  w.factory("cpi0_sig_y[1000,0,8000]")

  # now mis rec yield (ratio of this to sig should be the same for b0 and bs but will be different for g vs pi0)
  w.factory("misrec_to_sig_rat_g[0.2,0.001,0.6]")
  w.factory("misrec_to_sig_rat_pi0[0.2,0.001,0.6]")
  w.factory("prod::b0g_misrec_y( misrec_to_sig_rat_g, b0g_sig_y )")
  w.factory("prod::b0pi0_misrec_y( misrec_to_sig_rat_pi0, b0pi0_sig_y )")
  w.factory("prod::bsg_misrec_y( misrec_to_sig_rat_g, bsg_sig_y )")
  w.factory("prod::bspi0_misrec_y( misrec_to_sig_rat_pi0, bspi0_sig_y )")
  w.factory("prod::cg_misrec_y( misrec_to_sig_rat_g, cg_sig_y )")
  w.factory("prod::cpi0_misrec_y( misrec_to_sig_rat_pi0, cpi0_sig_y )")

  # the cases of B->D*pipi, B->D*KK, Lb->D*ph all involve a misID so will
  # be different for B0 and Bs (as they differ with a K or pi misID) however
  # for all of these the ratio of g -> pi0 should be the same
  # there is also Bs->D*KK which should scale the same for g and pi0 modes
  w.factory("misid_g2pi0_rat[0.1,0.0001,10.]")
  w.factory("b0g_bdstpp_y[1000,0,12000]")
  w.factory("bsg_bdstpp_y[1000,0,12000]")
  w.factory("prod::b0pi0_bdstpp_y( misid_g2pi0_rat, b0g_bdstpp_y )")
  w.factory("prod::bspi0_bdstpp_y( misid_g2pi0_rat, bsg_bdstpp_y )")
  w.factory("b0g_bdstkk_y[1000,0,12000]")
  w.factory("bsg_bdstkk_y[1000,0,12000]")
  w.factory("prod::b0pi0_bdstkk_y( misid_g2pi0_rat, b0g_bdstkk_y )")
  w.factory("prod::bspi0_bdstkk_y( misid_g2pi0_rat, bsg_bdstkk_y )")
  w.factory("b0g_lbdstph_y[1000,0,12000]")
  w.factory("bsg_lbdstph_y[1000,0,12000]")
  w.factory("prod::b0pi0_lbdstph_y( misid_g2pi0_rat, b0g_lbdstph_y )")
  w.factory("prod::bspi0_lbdstph_y( misid_g2pi0_rat, bsg_lbdstph_y )")
  w.factory("bsdstkk_to_bdstkk_rat[1.,0.1,2.]")
  w.factory("prod::b0g_bsdstkk_y( bsdstkk_to_bdstkk_rat, b0g_bdstkk_y )")
  w.factory("prod::b0pi0_bsdstkk_y( bsdstkk_to_bdstkk_rat, b0pi0_bdstkk_y )")
  w.factory("prod::bsg_bsdstkk_y( bsdstkk_to_bdstkk_rat, bsg_bdstkk_y )")
  w.factory("prod::bspi0_bsdstkk_y( bsdstkk_to_bdstkk_rat, bspi0_bdstkk_y )")

  # B -> DKpi same logic as misrec
  w.factory("bdkp_to_sig_rat_g[0.2,0.001,0.6]")
  w.factory("bdkp_to_sig_rat_pi0[0.2,0.001,0.6]")
  w.factory("prod::b0g_bdkp_y( bdkp_to_sig_rat_g, b0g_sig_y )")
  w.factory("prod::b0pi0_bdkp_y( bdkp_to_sig_rat_pi0, b0pi0_sig_y )")
  w.factory("prod::bsg_bdkp_y( bdkp_to_sig_rat_g, bsg_sig_y )")
  w.factory("prod::bspi0_bdkp_y( bdkp_to_sig_rat_pi0, bspi0_sig_y )")

  # infact can be much more sophisitcated with efficiencies etc here
  # i.e. if the PID cut distinguishing B -> D0 Kpi from B -> D0 pipi is binary one knows the exact yield of the cross feed in each
  # the B -> DKp cross feed yield (float for now)
  w.factory("cg_bdstkp_y[3000,0,12000]")
  w.factory("cpi0_bdstkp_y[800,0,4000]")

  # B -> D* K / B -> D* pi (adding random pi- to B0 and random K- to Bs0)
  # so ratio to signal should be same for both g and pi0 modes but
  # different for B0 -> Bs
  w.factory("bdsth_to_sig_rat_addpi[0.2,0.001,0.6]")
  w.factory("bdsth_to_sig_rat_addk[0.2,0.001,0.6]")
  w.factory("prod::b0g_bdsth_y( bdsth_to_sig_rat_addpi, b0g_sig_y )")
  w.factory("prod::b0pi0_bdsth_y( bdsth_to_sig_rat_addpi, b0pi0_sig_y )")
  w.factory("prod::bsg_bdsth_y( bdsth_to_sig_rat_addk, bsg_sig_y )")
  w.factory("prod::bspi0_bdsth_y( bdsth_to_sig_rat_addk, bspi0_sig_y )")
  # the B- -> D* pi- (can probably constrain this better as well)
  w.factory("cg_bdsth_y[3000,0,12000]")
  w.factory("cpi0_bdsth_y[800,0,4000]")

  # Lb -> Dph (mid-ID k for p and pi for p and add random g or pi0) so will be different for all 4 really
  # express this ratio to the Lb -> D*ph one (they should be similar in magnitude?)
  w.factory("lbdph_to_lbdstph_b0g[1.,0.5,2.]")
  w.factory("lbdph_to_lbdstph_b0pi0[1.,0.5,2.]")
  w.factory("lbdph_to_lbdstph_bsg[1.,0.5,2.]")
  w.factory("lbdph_to_lbdstph_bspi0[1.,0.5,2.]")
  w.factory("prod::b0g_lbdph_y( lbdph_to_lbdstph_b0g, b0g_lbdstph_y )")
  w.factory("prod::b0pi0_lbdph_y( lbdph_to_lbdstph_b0pi0, b0pi0_lbdstph_y )")
  w.factory("prod::bsg_lbdph_y( lbdph_to_lbdstph_bsg, bsg_lbdstph_y )")
  w.factory("prod::bspi0_lbdph_y( lbdph_to_lbdstph_bspi0, bspi0_lbdstph_y )")

  # Part reco shape should have same Bs / B0 ratio
  w.factory("partrec_to_sig_rat_g[0.2,0.001,0.6]")
  w.factory("partrec_to_sig_rat_pi0[0.2,0.001,0.6]")
  w.factory("prod::b0g_partrec_y( partrec_to_sig_rat_g, b0g_sig_y )")
  w.factory("prod::b0pi0_partrec_y( partrec_to_sig_rat_pi0, b0pi0_sig_y )")
  w.factory("prod::bsg_partrec_y( partrec_to_sig_rat_g, bsg_sig_y )")
  w.factory("prod::bspi0_partrec_y( partrec_to_sig_rat_pi0, bspi0_sig_y )")

  # make the yields (different for ctrl and b0 / bs)


  b_components = [ 'sig', 'misrec', 'bdstpp', 'bdstkk', 'bsdstkk', 'bdkp', 'bdsth', 'partrec', 'lbdph', 'lbdstph', 'comb' ]
  for samp in ['b0g','b0pi0','bsg','bspi0']:
    fact_str = "SUM::data_pdf_%s("%samp
    for comp in b_components:
      fact_str += "%s_%s_y*%s_mc_pdf_%s,"%(samp,comp,comp,samp)
    fact_str = fact_str[:-1] +")"
    w.factory(fact_str)
    w.pdf('data_pdf_%s'%samp).Print('v')

  ctrl_components = [ 'sig', 'misrec', 'bdstkp', 'bdsth', 'comb' ]
  for samp in ['cg','cpi0']:
    fact_str = "SUM::data_pdf_%s("%samp
    for comp in ctrl_components:
      fact_str += "%s_%s_y*%s_mc_pdf_%s,"%(samp,comp,comp,samp)
    fact_str = fact_str[:-1] +")"
    w.factory(fact_str)
    w.pdf('data_pdf_%s'%samp).Print('v')

  CreateSimPdf(w,'data')
  CreateSimData(w,'data')


  # Now fix appropriate parameters
  # To start with we'll fix all shape parameters from MC and just float the yields (and exponential slope)
  for comp in b_components:
    if comp == 'comb': continue # no pre-defined shape for combinatorial
    if comp == 'bsdstkk': continue # this params for this piece are covered by bdstkk
    w.set('%s_mc_sim_pdf_pars'%comp).setAttribAll("Constant")

  # Now relax the constraints on a few important params
  w.var("b0g_mean").setConstant(False)
  w.var("b0g_sigma").setConstant(False)
  #w.var("dm_b02bs").setConstant(False)
  #w.var("dm_g2pi0").setConstant(False)
  #w.var("ssig_b02bs").setConstant(False)
  #w.var("ssig_g2pi0").setConstant(False)

  #w.var("b0g_misrec_mean").setConstant(False)
  #w.var("b0g_misrec_sigma").setConstant(False)
  #w.var("dm_missg2addg").setConstant(False)
  #w.var("ssig_missg2addg").setConstant(False)

  w.Print('v')
  w.pdf('data_sim_pdf').Print('v')
  w.data('data_sim_data').Print('v')

  # free fit first
  if free:
    for i, samp in enumerate(samples):
      pdfname  = 'data_pdf_%s'%(samp)
      dsetname = 'data_%s'%(samp)
      w.pdf(pdfname).fitTo(w.data(dsetname)) # nothing to fit in this case
      pars = w.pdf(pdfname).getParameters(RooArgSet(w.var("B_DTFDict_D0_B_M")))
      w.saveSnapshot('data_free_fit_%s'%samp,pars)

  if sim:
    pdfname = 'data_sim_pdf'
    dsetname = 'data_sim_data'
    w.pdf(pdfname).fitTo( w.data(dsetname) )
    pars = w.pdf(pdfname).getParameters(RooArgSet(w.var("B_DTFDict_D0_B_M")))
    w.saveSnapshot('data_sim_fit',pars)

  w.writeToFile('files/w_ctrl.root')

def ControlDataPlot(free=True, sim=True):

  tf = TFile('files/w_ctrl.root')
  w = tf.Get('w_ctrl')

  samples = [ 'b0g', 'bsg', 'cg', 'b0pi0', 'bspi0', 'cpi0' ]

  cs = []

  # free fit
  if free:
    for i, samp in enumerate(samples):
      cs.append( TCanvas('c_data_free_%s'%samp,'',700,600) )
      pdfname  = 'data_pdf_%s'%(samp)
      dsetname = 'data_%s'%(samp)
      w.loadSnapshot('data_free_fit_%s'%samp)
      plotData( cs[-1], w, pdfname, dsetname, True, samp)
      cs[-1].Update()
      cs[-1].Modified()
      cs[-1].Print("plots/ctrl_data_free_fit_%s.pdf"%samp)
      plotData( cs[-1], w, pdfname, dsetname, True, samp, True)
      cs[-1].Update()
      cs[-1].Modified()
      cs[-1].Print("plots/ctrl_data_free_fit_%s_log.pdf"%samp)

  # sim fit
  if sim:
    w.loadSnapshot('data_sim_fit')
    for i, samp in enumerate(samples):
      cs.append( TCanvas('c_data_sim_%s'%samp,'',700,600) )
      pdfname  = 'data_pdf_%s'%(samp)
      dsetname = 'data_%s'%(samp)
      plotData( cs[-1], w, pdfname, dsetname, True, samp)
      cs[-1].Update()
      cs[-1].Modified()
      cs[-1].Print("plots/ctrl_data_sim_fit_%s.pdf"%samp)
      plotData( cs[-1], w, pdfname, dsetname, True, samp, True)
      cs[-1].Update()
      cs[-1].Modified()
      cs[-1].Print("plots/ctrl_data_sim_fit_%s_log.pdf"%samp)

