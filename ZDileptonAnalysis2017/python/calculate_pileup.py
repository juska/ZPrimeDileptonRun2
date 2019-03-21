#chad harrington 2019
#execute: python python/calculate_pileup.py -mc dir1:file1 dir2:file2 -d dir:file
#default: python python/calculate_pileup.py
#example: python python/calculate_pileup.py -mc root://131.225.204.161:1094//store/user/cihar29/2017Analysis_trees:TTjets root://131.225.204.161:1094//store/user/cihar29/2017Analysis_trees:ww -d /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/root_trees:mu_Data
import argparse
from files import mcTups
import ROOT
ROOT.gROOT.SetBatch(True)

def mcList_default() :
  mcList = []
  for mcTup in mcTups : mcList += [(mcTup.dir, mcTup.name)]
  return mcList

def data_default() :
  return [("/uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017","mu_Data")]

def arg2List( args ) :
  list = []
  for pair in args :
    split = pair.split(":")
    l = len(split)
    if l < 2 : raise Exception("Must provide a directory and dataset, as in dir:dataset")

    dir, name = ":".join(split[:l-1]), split[l-1] #get last ':'
    #print dir, name
    list += [(dir, name)]
  return list

def main() :
  parser = argparse.ArgumentParser()
  parser.add_argument( "-mc", "--mcList", nargs='+', required=False,
                       help="dir1:mcDataset1 dir2:mcDataset2"
                     )
  parser.add_argument( "-d", "--data", action='append', required=False,
                       help="dir:data where dataNOMINAL, dataUP, and dataDOWN files exist"
                     )
  args = parser.parse_args()

  if args.mcList is None : mc = mcList_default()
  else :                   mc = arg2List(args.mcList)
  if args.data is None : data = data_default()
  else :                 data = arg2List(args.data)

  calc_pileup( mc, data )

def calc_pileup( mcList=mcList_default(), data=data_default() ) :

  ddir, dname = data[0][0], data[0][1]
  wgts = ("NOMINAL", "UP", "DOWN")
  dList = []
  for wgt in wgts : dList += [(ddir, dname+wgt)]

  dDict  = getHists( dList )
  mcDict = getHists( mcList )

  drawPileup( dDict[dname+"NOMINAL"], mcDict )

  drawWeights( dname, wgts, dDict, mcDict )

def drawPileup( h_data, mcDict ) :
  setStyle()
  c = ROOT.TCanvas("c", "c", 600, 600)

  leg = getLegend( .5, .9-(len(mcDict)+1)*0.03, .7, .9 )

  h_data.GetXaxis().SetTitle("#mu")
  h_data.SetMarkerSize(0.7)
  h_data.SetMarkerStyle(4)
  h_data.Draw()
  leg.AddEntry(h_data, "%-15s (Mean %4.1f, RMS %4.1f)" % ("data", h_data.GetMean(), h_data.GetRMS()), "P")

  for name, hist in mcDict.iteritems() :
    leg.AddEntry(hist, "%-15s (Mean %4.1f, RMS %4.1f)" % (name, hist.GetMean(), hist.GetRMS()), "L")
    hist.Draw("histsame")

  leg.Draw()
  drawText()
  c.Print("mu.pdf")

def drawWeights( dname, wgts, dDict, mcDict ) :
  setStyle()
  c = ROOT.TCanvas("c", "c", 600, 600)

  leg = getLegend( .6, .9-(3*len(mcDict))*0.03, .8, .9 )

  h_data = dDict[dname+"NOMINAL"].Clone("h_data_wgt")
  h_data.Divide( h_data )
  h_data.GetYaxis().SetRangeUser(0, 5)
  h_data.Draw("hist p")

  outName = "mu_weights.root"
  outFile = ROOT.TFile(outName, "RECREATE")
  outFile.cd()

  wgthists = []
  istyle = 1
  for wgt in wgts :
    outFile.mkdir( wgt )
    outFile.cd( outName + ":/" + wgt )

    icolor = 1
    for name, hist in mcDict.iteritems() :
      wgtname = name + "_" + wgt

      wgthists += [dDict[ dname+wgt ].Clone(wgtname)]
      wgthists[-1].Divide( hist )

      wgthists[-1].SetLineColor(icolor)
      icolor += 1
      wgthists[-1].SetLineStyle(istyle)
      wgthists[-1].Draw("histsame")

      leg.AddEntry(wgthists[-1], wgtname, "L")
      wgthists[-1].Write(wgtname)
    istyle += 3
  outFile.Close()

  leg.Draw()
  drawText( "Weights" )
  c.Print("mu_weights.pdf")

def getHists( list ) :
  icolor = 1
  hists = {}
  for dir, name in list :

    file = ROOT.TFile.Open( "%s/%s.root" % (dir, name) )
    if file == None : exit()
    ROOT.TH1.AddDirectory(0)
    hists[name] = file.Get( "fullMu" )
    hists[name].Scale( 1 / hists[name].Integral() )

    hists[name].SetLineColor( icolor )
    icolor += 1
  return hists

def drawText( rightText="Run 2017 - 41.9 fb^{-1} (13 TeV)" ) :

  text = ROOT.TLatex()
  text.SetNDC()

  text.SetTextSize(0.035)
  text.SetTextFont(42)
  l = len(rightText)
  x = 1-l/68.
  if len(rightText) < 10 : x = 1-l/40.
  text.DrawLatex(x, 0.96, rightText)

  text.SetTextSize(0.05)
  text.SetTextFont(61)
  text.DrawLatex(0.12, 0.96, "CMS")
  text.SetTextSize(0.03)

def getLegend( x1, y1, x2, y2 ) :
  leg = ROOT.TLegend( x1, y1, x2, y2 )
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.025)
  leg.SetTextFont(42)
  return leg

def setStyle() :

  ROOT.gStyle.SetPadTopMargin(0.05)
  ROOT.gStyle.SetPadBottomMargin(0.1)
  ROOT.gStyle.SetPadLeftMargin(0.1)
  ROOT.gStyle.SetPadRightMargin(0.04)

  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)

  ROOT.gStyle.SetTitleFont(42, "XYZ")
  ROOT.gStyle.SetTitleSize(0.06, "XYZ")
  ROOT.gStyle.SetTitleXOffset(1.1)
  ROOT.gStyle.SetTitleYOffset(1.4)

  ROOT.gStyle.SetLabelFont(42, "XYZ")
  ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
  ROOT.gStyle.SetLabelSize(0.04, "XYZ")

  ROOT.gStyle.SetPadTickX(1)
  ROOT.gStyle.SetPadTickY(1)

if __name__ == "__main__":
    main()
