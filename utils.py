from ROOT import gROOT, TPad, TLine, TH1F, TLatex, TLegend
from ROOT import kBlack, kBlue, kRed, kGreen, kOrange, kMagenta, kPink, kYellow, kCyan
from ROOT import RooFit as rf
from ROOT import RooLinkedList, RooSimultaneous, RooDataSet, RooArgSet
import roottex as rt

# Shared params to describe mass shift and width scale from g -> pi0
def SharedParams(w):
  w.factory("dm_g2pi0[5,-5.,20]")
  w.factory("ssig_g2pi0[1.,0.5,2.]")
  w.factory("dm_b02bs[70,50,100]")
  w.factory("ssig_b02bs[1.,0.1,2.]")

# Make sim pdf
def CreateSimPdf(w, typ):
  sim_pdf = RooSimultaneous("%s_sim_pdf"%typ,"",w.cat("cat"))
  # loop on cats
  found_one = False
  for c in range(w.cat("cat").numTypes()):
    w.cat("cat").setIndex(c)
    samp = w.cat("cat").getLabel()
    subpdfname = '%s_pdf_%s'%(typ,samp)
    if w.pdf(subpdfname):
      sim_pdf.addPdf( w.pdf(subpdfname), samp )
      found_one = True
  if found_one:
    getattr(w,'import')(sim_pdf)

# Make sim data (this needs a bit of hacking because of pyROOT)
def CreateSimData(w, typ):
  sim_dsets = []
  for c in range(w.cat("cat").numTypes()):
    w.cat("cat").setIndex(c)
    samp = w.cat("cat").getLabel()
    name = "%s_sim_data"%typ if c==0 else "%s_sim_data%d"%(typ,c)
    sim_dsets.append( RooDataSet(name,"",RooArgSet(w.var("B_DTFDict_D0_B_M")),rf.Index(w.cat("cat")),rf.Import(samp,w.data("%s_%s"%(typ,samp)))) )
    if c!=0:
      sim_dsets[0].append(sim_dsets[c])
  getattr(w,'import')(sim_dsets[0])

def plot(pad, w, pdfname, dsetname, resid=False, cat=None, popts=None, dopts=None):

  if cat:
    if cat.startswith('b0_'): cat = cat.replace('b0_','')
    if cat.startswith('bs_'): cat = cat.replace('bs_','')

  xtitle = { 'b0g'   : '#it{m}('+rt.Dzb+rt.g+rt.Kp+rt.pim  +') [MeV/c^{2}]' ,
             'b0pi0' : '#it{m}('+rt.Dzb+rt.piz+rt.Kp+rt.pim+') [MeV/c^{2}]',
             'bsg'   : '#it{m}('+rt.Dzb+rt.g+rt.Km+rt.pip  +') [MeV/c^{2}]' ,
             'bspi0' : '#it{m}('+rt.Dzb+rt.piz+rt.Km+rt.pip+') [MeV/c^{2}]'
           }
  leg    = { 'b0g'   : '#it{B}^{0} #rightarrow #bar{#it{D}}^{0}#it{#gamma}#it{K}^{+}#it{#pi}^{-}',
             'b0pi0' : '#it{B}^{0} #rightarrow #bar{#it{D}}^{0}#it{#pi}^{0}#it{K}^{+}#it{#pi}^{-}',
             'bsg'   : '#it{B_{s}}^{#kern[-0.2]{0}} #rightarrow #bar{#it{D}}^{0}#it{#gamma}#it{K}^{-}#it{#pi}^{+}',
             'bspi0' : '#it{B_{s}}^{#kern[-0.2]{0}} #rightarrow #bar{#it{D}}^{0}#it{#pi}^{0}#it{K}^{-}#it{#pi}^{+}'
           }

  if not popts: popts = RooLinkedList()
  if not dopts: dopts = RooLinkedList()

  pl = w.var('B_DTFDict_D0_B_M').frame()
  if cat: pl.GetXaxis().SetTitle( xtitle[cat] )
  if w.data(dsetname): w.data(dsetname).plotOn(pl,dopts)
  if w.pdf(pdfname): w.pdf(pdfname).plotOn(pl,popts)
  if resid and w.data(dsetname) and w.pdf(pdfname):
    pu = TPad(pad.GetName()+'_pu','',0., 0.33, 1., 1.)
    pd = TPad(pad.GetName()+'_pd','',0., 0., 1., 0.33)
    pu.SetBottomMargin(0.03)
    pd.SetTopMargin(0.03)
    pd.SetBottomMargin(0.35)
    pad.cd()
    pu.Draw()
    pd.Draw()
    pu.cd()
    pl.GetYaxis().SetTitleOffset(1.1)
    pl.GetXaxis().SetLabelSize(0.)
    pl.Draw()
    pd.cd()
    pull = pl.pullHist()
    # roofit does something funny with pull hist range so need a dummy hist for it
    pull.Draw()
    ymax = max( abs(pull.GetYaxis().GetXmin()), pull.GetYaxis().GetXmax() )
    pulld = TH1F( pull.GetName()+'_pd','', pl.GetNbinsX(), pl.GetXaxis().GetXmin(), pl.GetXaxis().GetXmax() )
    pulld.GetYaxis().SetRangeUser(-ymax,ymax)
    pulld.GetXaxis().SetTitle( pl.GetXaxis().GetTitle() )
    pulld.GetYaxis().SetTitle("Pull")
    pl.GetXaxis().SetTitle("")
    pulld.GetXaxis().SetTitleSize(0.13)
    pulld.GetYaxis().SetTitleSize(0.13)
    pulld.GetXaxis().SetLabelSize(0.13)
    pulld.GetYaxis().SetLabelSize(0.13)
    pulld.GetXaxis().SetTitleOffset(1.1)
    pulld.GetYaxis().SetTitleOffset(0.6)
    pulld.GetYaxis().SetNdivisions(505)
    pulld.GetXaxis().SetTickLength(0.06)
    ll = TLine()
    ll.SetLineWidth(4)
    ll.SetLineColor(4)
    pd.cd()
    pd.Clear()
    pulld.Draw("HIST")
    pull.Draw("EPsame")
    ll.DrawLine( pulld.GetXaxis().GetXmin(), 0., pulld.GetXaxis().GetXmax(), 0. )
    pulld.Draw("AXISsame")
    gROOT.Append(pulld)
  else:
    pad.cd()
    pl.Draw()

  pad.cd()
  if cat:
    lat = TLatex()
    lat.SetTextSize(0.05)
    lat.SetNDC()
    lat.DrawLatex(0.6,0.8,leg[cat])
    gROOT.Append(lat)
  pad.Update()
  pad.Modified()

def plotData(pad, w, pdfname, dsetname, resid=False, cat=None, log=False):

  if cat.startswith('b0_'): cat = cat.replace('b0_','')
  if cat.startswith('bs_'): cat = cat.replace('bs_','')

  xtitle = { 'b0g'   : '#it{m}('+rt.Dzb+rt.g+rt.Kp+rt.pim  +') [MeV/c^{2}]' ,
             'b0pi0' : '#it{m}('+rt.Dzb+rt.piz+rt.Kp+rt.pim+') [MeV/c^{2}]' ,
             'bsg'   : '#it{m}('+rt.Dzb+rt.g+rt.Km+rt.pip  +') [MeV/c^{2}]' ,
             'bspi0' : '#it{m}('+rt.Dzb+rt.piz+rt.Km+rt.pip+') [MeV/c^{2}]' ,
             'cg'    : '#it{m}('+rt.Dzb+rt.g+rt.pip+rt.pip  +') [MeV/c^{2}]',
             'cpi0'  : '#it{m}('+rt.Dzb+rt.piz+rt.pim+rt.pip+') [MeV/c^{2}]'
           }
  leg    = { 'b0g'   : { 'sig'    : rt.Decay('Bd','Dzb','g','Kp','pim'),
                         'misrec' : rt.Decay('Bd','Dzb','piz','Kp','pim'),
                         'bdstpp' : rt.Decay('Bd','Dzb','g','pip','pim'),
                         'bdstkk' : rt.Decay('Bd','Dzb','g','Kp','Km'),
                         'bsdstkk': rt.Decay('Bs','Dzb','g','Kp','Km'),
                         'bdkp'   : rt.Decay('Bd','Dzb','Kp','pim'),
                         'bdsth'  : rt.Decay('Bm','Dzb','g','pim'),
                         'partrec': 'Part Reco',
                         'lbdph'  : rt.Decay('Lb','Dzb','p','pim'),
                         'lbdstph': rt.Decay('Lb','Dzb','g','p','pim'),
                         'comb'   : 'Combinatorial'
                       },
             'b0pi0' : { 'sig'    : rt.Decay('Bd','Dzb','piz','Kp','pim'),
                         'misrec' : rt.Decay('Bd','Dzb','g','Kp','pim'),
                         'bdstpp' : rt.Decay('Bd','Dzb','piz','pip','pim'),
                         'bdstkk' : rt.Decay('Bd','Dzb','piz','Kp','Km'),
                         'bsdstkk': rt.Decay('Bs','Dzb','piz','Kp','Km'),
                         'bdkp'   : rt.Decay('Bd','Dzb','Kp','pim'),
                         'bdsth'  : rt.Decay('Bm','Dzb','piz','pim'),
                         'partrec': 'Part Reco',
                         'lbdph'  : rt.Decay('Lb','Dzb','p','pim'),
                         'lbdstph': rt.Decay('Lb','Dzb','piz','p','pim'),
                         'comb'   : 'Combinatorial'
                        },
             'bsg'   : { 'sig'    : rt.Decay('Bs','Dzb','g','Km','pip'),
                         'misrec' : rt.Decay('Bs','Dzb','piz','Km','pip'),
                         'bdstpp' : rt.Decay('Bs','Dzb','g','pim','pip'),
                         'bdstkk' : rt.Decay('Bs','Dzb','g','Km','Kp'),
                         'bsdstkk': rt.Decay('Bs','Dzb','g','Km','Kp'),
                         'bdkp'   : rt.Decay('Bs','Dzb','Km','pip'),
                         'bdsth'  : rt.Decay('Bp','Dzb','g','pip'),
                         'partrec': 'Part Reco',
                         'lbdph'  : rt.Decay('Lb','Dzb','p','Km'),
                         'lbdstph': rt.Decay('Lb','Dzb','g','p','Km'),
                         'comb'   : 'Combinatorial'
                        },
             'bspi0' : { 'sig'    : rt.Decay('Bs','Dzb','piz','Km','pip'),
                         'misrec' : rt.Decay('Bs','Dzb','g','Km','pip'),
                         'bdstpp' : rt.Decay('Bs','Dzb','piz','pim','pip'),
                         'bdstkk' : rt.Decay('Bs','Dzb','piz','Km','Kp'),
                         'bsdstkk': rt.Decay('Bs','Dzb','piz','Km','Kp'),
                         'bdkp'   : rt.Decay('Bs','Dzb','Km','pip'),
                         'bdsth'  : rt.Decay('Bp','Dzb','piz','pip'),
                         'partrec': 'Part Reco',
                         'lbdph'  : rt.Decay('Lb','Dzb','p','Km'),
                         'lbdstph': rt.Decay('Lb','Dzb','piz','p','Km'),
                         'comb'   : 'Combinatorial'
                       },
             'cg'    : { 'sig'    : rt.Decay('Bd','Dzb','g','pip','pim'),
                         'misrec' : rt.Decay('Bd','Dzb','piz','pip','pim'),
                         'bdstkp' : rt.Decay('Bs','Dzb','g','Km','pip'),
                         'bdkp'   : rt.Decay('Bd','Dzb','Kp','pim'),
                         'bdsth'  :  rt.Decay('Bm','Dzb','g','pim'),
                         'comb'   : 'Combinatorial'
                        },
              'cpi0' : { 'sig'    : rt.Decay('Bd','Dzb','piz','pip','pim'),
                         'misrec' : rt.Decay('Bd','Dzb','g','pip','pim'),
                         'bdstkp' : rt.Decay('Bs','Dzb','piz','Km','pip'),
                         'bdkp'   : rt.Decay('Bd','Dzb','Kp','pim'),
                         'bdsth'  : rt.Decay('Bm','Dzb','piz','pim'),
                         'comb'   : 'Combinatorial'
                       },
             }

  pl = w.var('B_DTFDict_D0_B_M').frame()
  if cat: pl.GetXaxis().SetTitle( xtitle[cat] )
  tleg = TLegend(0.6,0.5,0.9,0.95)
  tleg.SetFillColorAlpha(0,0.)
  tleg.SetLineColor(0)
  if w.data(dsetname):
    w.data(dsetname).plotOn(pl)
  if w.pdf(pdfname):
    w.pdf(pdfname).plotOn(pl)
    w.pdf(pdfname).plotOn(pl, rf.Components('comb_mc_pdf_%s'%cat)   , rf.LineStyle(2), rf.LineColor(kBlack) )
    tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['comb'], "L" )
    if cat!='cg' and cat !='cpi0':
      w.pdf(pdfname).plotOn(pl, rf.Components('partrec_mc_pdf_%s'%cat), rf.LineStyle(2), rf.LineColor(kRed+3) )
      tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['partrec'], "L" )
      w.pdf(pdfname).plotOn(pl, rf.Components('lbdph_mc_pdf_%s'%cat)  , rf.LineStyle(2), rf.LineColor(kMagenta) )
      tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['lbdph'], "L" )
      w.pdf(pdfname).plotOn(pl, rf.Components('lbdstph_mc_pdf_%s'%cat), rf.LineStyle(2), rf.LineColor(kMagenta+3) )
      tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['lbdstph'], "L" )
    w.pdf(pdfname).plotOn(pl, rf.Components('bdsth_mc_pdf_%s'%cat)  , rf.LineStyle(2), rf.LineColor(kRed-3) )
    tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['bdsth'], "L" )
    if cat!='cg' and cat !='cpi0':
      w.pdf(pdfname).plotOn(pl, rf.Components('bdkp_mc_pdf_%s'%cat)   , rf.LineStyle(2), rf.LineColor(kOrange-3) )
      tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['bdkp'], "L" )
      w.pdf(pdfname).plotOn(pl, rf.Components('bsdstkk_mc_pdf_%s'%cat), rf.LineStyle(2), rf.LineColor(kGreen+1) )
      tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['bsdstkk'], "L" )
      w.pdf(pdfname).plotOn(pl, rf.Components('bdstkk_mc_pdf_%s'%cat) , rf.LineStyle(2), rf.LineColor(kGreen+3) )
      tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['bdstkk'], "L" )
      w.pdf(pdfname).plotOn(pl, rf.Components('bdstpp_mc_pdf_%s'%cat) , rf.LineStyle(2), rf.LineColor(kRed) )
      tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['bdstpp'], "L" )
    else:
      w.pdf(pdfname).plotOn(pl, rf.Components('bdstkp_mc_pdf_%s'%cat)   , rf.LineStyle(2), rf.LineColor(kOrange-3) )
      tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['bdstkp'], "L" )
    w.pdf(pdfname).plotOn(pl, rf.Components('misrec_mc_pdf_%s'%cat) , rf.LineStyle(2), rf.LineColor(kCyan) )
    tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['misrec'], "L" )
    w.pdf(pdfname).plotOn(pl, rf.Components('sig_mc_pdf_%s'%cat)    , rf.LineStyle(2), rf.LineColor(kBlue+3) )
    tleg.AddEntry( pl.getObject( int(pl.numItems())-1 ), leg[cat]['sig'], "L" )
    w.pdf(pdfname).plotOn(pl)
  if resid and w.data(dsetname) and w.pdf(pdfname):
    pu = TPad(pad.GetName()+'_pu','',0., 0.33, 1., 1.)
    pd = TPad(pad.GetName()+'_pd','',0., 0., 1., 0.33)
    pu.SetBottomMargin(0.03)
    pd.SetTopMargin(0.03)
    pd.SetBottomMargin(0.35)
    pad.cd()
    pu.Draw()
    pd.Draw()
    pu.cd()
    pl.GetYaxis().SetTitleOffset(1.1)
    pl.GetXaxis().SetLabelSize(0.)
    pl.Draw()
    if log:
      pl.GetYaxis().SetRangeUser(0.1, 2.*pl.GetMaximum())
      pl.Draw()
      pu.SetLogy()
    pd.cd()
    pull = pl.pullHist()
    # roofit does something funny with pull hist range so need a dummy hist for it
    pull.Draw()
    ymax = max( abs(pull.GetYaxis().GetXmin()), pull.GetYaxis().GetXmax() )
    pulld = TH1F( pull.GetName()+'_pd','', pl.GetNbinsX(), pl.GetXaxis().GetXmin(), pl.GetXaxis().GetXmax() )
    pulld.GetYaxis().SetRangeUser(-ymax,ymax)
    pulld.GetXaxis().SetTitle( pl.GetXaxis().GetTitle() )
    pulld.GetYaxis().SetTitle("Pull")
    pl.GetXaxis().SetTitle("")
    pulld.GetXaxis().SetTitleSize(0.13)
    pulld.GetYaxis().SetTitleSize(0.13)
    pulld.GetXaxis().SetLabelSize(0.13)
    pulld.GetYaxis().SetLabelSize(0.13)
    pulld.GetXaxis().SetTitleOffset(1.1)
    pulld.GetYaxis().SetTitleOffset(0.6)
    pulld.GetYaxis().SetNdivisions(505)
    pulld.GetXaxis().SetTickLength(0.06)
    ll = TLine()
    ll.SetLineWidth(4)
    ll.SetLineColor(4)
    pd.cd()
    pd.Clear()
    pulld.Draw("HIST")
    pull.Draw("EPsame")
    ll.DrawLine( pulld.GetXaxis().GetXmin(), 0., pulld.GetXaxis().GetXmax(), 0. )
    pulld.Draw("AXISsame")
    gROOT.Append(pulld)
  else:
    pad.cd()
    pl.Draw()
    if log:
      pl.GetYaxis().SetRangeUser(0.1, 2*pl.GetMaximum())
      pl.Draw()
      pad.SetLogy()

  gROOT.Append(tleg)
  pad.cd()
  if not log: tleg.Draw("same")
  pad.Update()
  pad.Modified()
