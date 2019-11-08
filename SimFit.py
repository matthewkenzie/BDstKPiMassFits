from argparse import ArgumentParser
parser = ArgumentParser(description="B2DstKpi Mass Fits")
parser.add_argument("-m","--mc"  ,default=False,action="store_true",help='Rerun the MC shape fits')
parser.add_argument("-d","--data",default=False,action="store_true",help='Run the Data fit (no control channels)')
parser.add_argument("-c","--ctrl",default=False,action="store_true",help='Run the data fit with control channels')
parser.add_argument("-a","--all", default=False,action="store_true",help='Run the full chain')
args = parser.parse_args()
if args.all:
  args.mc = True
  args.data = True
  args.ctrl = True

import sys
import os
os.system('mkdir -p plots')
os.system('mkdir -p files')

from ROOT import gROOT, gStyle
gROOT.SetBatch()
gROOT.ProcessLine(".x lhcbStyle.C")
gStyle.SetTitleXSize(-1.07);
gStyle.SetTitleYSize(0.07);
gStyle.SetTitleXOffset(0.9);
gStyle.SetTitleYOffset(0.9);
gStyle.SetPadLeftMargin(0.15);

from mc_shapes import *

if args.mc:

  ### Signal Shapes ###
  SignalShape()
  ### MisRec Shapes ###
  MisRecShape()
  ### MisID Shapes
  BDstPPShape()
  BDstKKShape()
  ### Part Comb Shapes ###
  BDKPShape(False)
  BDstHShape(False)
  ### Part Reco Shape ###
  PartRecShape()
  ### Lb Shapes ###
  Lb2DPHShape()
  Lb2DstPHShape()

from data_fit import DataFit, DataPlot
if args.data:
  ### Data Fit ###
  DataFit()
  DataPlot()

from data_fit import ControlDataFit, ControlDataPlot
if args.ctrl:
  ### Control Data Fit ###
  ControlDataFit()
  ControlDataPlot()
