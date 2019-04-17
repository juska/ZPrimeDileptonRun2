# chad harrington and bahareh roozbahani 2019
# python -u python/plot_all.py -drs ON OFF -c mm em -cuts 4 7 > & log.txt &

import os
import argparse
import time
from collections import namedtuple

from files import mcTups, dataDict

drs_default = ["OFF", "ON", "ONbt", "ONnb", "REVERSE"]
channels_default = ["mm", "em", "ee"]
cuts_default = ["0", "1", "2", "3", "4", "5", "6", "7"]

def handleArgs() :
  parser = argparse.ArgumentParser()
  parser.add_argument( "-drs", "--drCuts", nargs='+', required=False,
                       help="drCuts to run over, like -drs OFF ON",
                       default=drs_default
                     )
  parser.add_argument( "-c", "--channels", nargs='+', required=False,
                       help="channels to run over, like -c mm",
                       default=channels_default
                     )
  parser.add_argument( "-cuts", "--cuts", nargs='+', required=False,
                       help="cuts to plot, like -cuts 7",
                       default=cuts_default
                     )
  parser.add_argument( "-o", "--outdir", type=str, action='store', required=False,
                       help="output plot directory",
                       default="plots"
                     )
  return parser.parse_args()

def main() :
  start = time.time()
  args = handleArgs()

  drs = args.drCuts
  channels = args.channels
  cuts = args.cuts

  print "drCuts:", " ".join(drs)
  print "channels:", " ".join(channels)
  print "cuts:", " ".join(cuts)

  channel_dict = {
                  'mm' : ('muon',                  '/uscms_data/d3/broozbah/ZPRIME_2017/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/OFF'),
                  'ee' : ('DoubleEG_2017_31Mar18', '/uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/OFF'),
                  'em' : ('muon',                  '/uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/OFF')
                 }

  plotTuple = namedtuple('plotTuple', 'hname xtitle xmin xmax')

  plotTups = (
    plotTuple( 'nbtag',       'Number of btagged Jets',                           0,  5    ),
    plotTuple( 'jet0btag',    'btag_{Leading Jet}',                               0,  1    ),
    plotTuple( 'jet1btag',    'btag_{Subleading Jet}',                            0,  1    ),
    plotTuple( 'jet0eta',     '#eta_{Leading Jet}',                               -5, 5    ),
    plotTuple( 'jet1eta',     '#eta_{Subleading Jet}',                            -5, 5    ),
    plotTuple( 'jet0pt',      'Leading Jet p_{T} [GeV]',                          0,  2000 ),
    plotTuple( 'jet1pt',      'Subleading Jet p_{T} [GeV]',                       0,  2000 ),
    plotTuple( 'minjet0pt',   'Jet p_{T}^{rmin leading lepton} [GeV]',            0,  1000 ),
    plotTuple( 'minjet1pt',   'Jet p_{T}^{rmin subleading lepton} [GeV]',         0,  1000 ),
    plotTuple( 'cleanjet0pt', 'Jet p_{T}^{cleaned from leading lepton} [GeV]',    0,  1000 ),
    plotTuple( 'cleanjet1pt', 'Jet p_{T}^{cleaned from subleading lepton} [GeV]', 0,  1000 ),
    plotTuple( 'jethT',       'H_{T} [GeV]',                                      0,  3000 ),
    plotTuple( 'nGoodJet',    'N_{Good Jets}',                                    0,  25   ),

    plotTuple( 'lep0eta',     '#eta_{Leading Lepton}',                  -5, 5    ),
    plotTuple( 'lep1eta',     '#eta_{Subleading Lepton}',               -5, 5    ),
    plotTuple( 'lep0pt',      'Leading Lepton p_{T} [GeV]',             0,  1000 ),
    plotTuple( 'lep1pt',      'Subleading Lepton p_{T} [GeV]',          0,  600  ),
    plotTuple( 'lep0perp',    'Leading Lepton p_{T}^{rel} [GeV]',       0,  200  ),
    plotTuple( 'lep1perp',    'Subleading Lepton p_{T}^{rel} [GeV]',    0,  200  ),
    plotTuple( 'lep0perp_in', 'Leading Lepton p_{T}^{rel} In [GeV]',    0,  200  ),
    plotTuple( 'lep1perp_in', 'Subleading Lepton p_{T}^{rel} In [GeV]', 0,  200  ),
    plotTuple( 'lepept',      'electron p_{T} [GeV]',                   0,  1000 ),
    plotTuple( 'lepmpt',      'muon p_{T} [GeV]',                       0,  1000 ),
    plotTuple( 'nGoodEle',    'N_{Good Electrons}',                     0,  5    ),
    plotTuple( 'nGoodMuon',   'N_{Good Muons}',                         0,  5    ),
    plotTuple( 'nPV',         'N_{Good Primary vertices}',              0,  80   ),

    plotTuple( 'metpt',     'MET p_{T}^{uncorr} [GeV]', 0, 1000 ),
    plotTuple( 'metcorrpt', 'MET p_{T} [GeV]',          0, 1000 ),
    plotTuple( 'dilepmass', 'M_{ll} [GeV]',             0, 1000 ),
    plotTuple( 'sT',        'H_{T}^{L} [GeV]',          0, 5000 ),
    plotTuple( 'sT_met',    'S_{T} [GeV]',              0, 5000 ),
    plotTuple( 'masslljjm', 'M_{lljjmet} [GeV]',        0, 5000 ),
    plotTuple( 'masslmin0', 'M_{leading lep,rmin jet} [GeV]',    0, 2000 ),
    plotTuple( 'masslmin1', 'M_{subleading lep,rmin jet} [GeV]', 0, 2000 ),

    plotTuple( 'rl0l1',        '#DeltaR(leading lepton, subleading lepton)', 0, 5  ),
    plotTuple( 'rl0cleanj',    '#DeltaR(leading lepton, cleaned jet)',       0, 5  ),
    plotTuple( 'rl1cleanj',    '#DeltaR(subleading lepton, cleaned jet)',    0, 5  ),
    plotTuple( 'rmin0',        '#DeltaR_{min}(leading lepton, jet)',         0, 5  ),
    plotTuple( 'rmin1',        '#DeltaR_{min}(subleading lepton, jet)',      0, 5  ),
    plotTuple( 'sumrmin',      '#DeltaR_{min0} + #DeltaR_{min1}',            0, 10 ),
    plotTuple( 'rbl',          '#DeltaR(b quark, lepton)',                   0, 5  ),
    plotTuple( 'dphi_jet0met', '#Delta #phi_{Leading Jet,MET}',              0, 3  ),
    plotTuple( 'dphi_jet1met', '#Delta #phi_{Subleading Jet,MET}',           0, 3  ),

    plotTuple( 'MT2r',           'MT2r',           0, 400 ),
    plotTuple( 'T_xl',           'T_xl',           0, 2   ),
    plotTuple( 'ANTIT_xl',       'ANTIT_xl',       0, 2   ),
    plotTuple( 'cosTheta1r',     'cosTheta1r',     -1, 1  ),
    plotTuple( 'cosTheta2r',     'cosTheta2r',     -1, 1  ),
    plotTuple( 'cosTheta1r_mt2', 'cosTheta1r_mt2', -1, 1  ),
    plotTuple( 'cosTheta2r_mt2', 'cosTheta2r_mt2', -1, 1  ),

    plotTuple( 'MT2s',       'MT2s',       0, 400 ),
    plotTuple( 'TOP_xl',     'TOP_xl',     0, 2   ),
    plotTuple( 'ANTITOP_xl', 'ANTITOP_xl', 0, 2   ),
    plotTuple( 'cosTheta1',  'cosTheta1',  -1, 1  ),
    plotTuple( 'cosTheta2',  'cosTheta2',  -1, 1  )
             )

  subplot = 'ratio'
  if subplot == 'ratio' : subymin = '0'
  else :                  subymin = '-2'

  os.system( "mkdir " + args.outdir )
  for drcut in drs :
    if drcut not in drs_default : print drcut + " not found!"; continue
    os.system( "mkdir " + args.outdir + "/" + drcut )

    for channel in channels :
      if channel not in channels_default : print channel + " not found!"; continue
      os.system( "mkdir " + args.outdir + "/" + drcut + "/" + channel )

      for cutNum in cuts :
        if cutNum not in cuts_default : print cutNum + " not found!"; continue
        os.system( "mkdir " + args.outdir + "/" + drcut + "/" + channel + "/" + cutNum )

        for plotTup in plotTups :

          logy = 'false'
          if 'pt' in plotTup.hname or 'sT' in plotTup.hname or plotTup.hname == 'jethT' or plotTup.hname == 'masslljjm' :
            logy = 'true'

          rebin = '2'
          if 'Good' in plotTup.hname or 'perp' in plotTup.hname or plotTup.hname == 'nbtag' :
            rebin = '1'

          writePars( channel_dict[channel][0], channel_dict[channel][1], args.outdir, channel, cutNum,
                     plotTup.hname, logy, plotTup.xmin, plotTup.xmax, plotTup.xtitle, subplot, subymin, rebin
                   )
          os.system( "plot plot_pars.txt" )

def writePars( dfile, dir, odir, channel, cut, hname, logy, xmin, xmax, xtitle, subplot, subymin, rebin ) :

  lines=[ "dataFileName   {}".format( dfile ),
          "plotData       true",
          "dir            {}/{}/".format( dir, channel ),
          "channel        {}".format( channel ),
          "fileNames      TTbar0-700 TTbar700-1000_total TTbar1000-inf_total TTsemilep wjets ww wz zz DYJetsToLL_M10to50 DYJetsToLL_M50 ST_t-channel_top_4f ST_t-channel_antitop_4f ST_tW_top_5f_inclusiveDecays_v2 ST_tW_antitop_5f_v2",
          "zprime         ZprimeToTT_M3000_W300_total",
          "gluon          RSGluonToTT_M-3000",
          "sigScale       10",
          "leftText       CMS",
          "rightText      Run 2017 - 41.9 fb^{-1} (13 TeV)",
          "logx           false",
          "logy           {}".format( logy ),
          "subplot        {}".format( subplot ),
          "outName        {}/{}/{}/{}".format( odir, channel, cut, hname ),
          "hname          {}_{}".format( cut, hname ),
          "xmin           {}".format( xmin ),
          "xmax           {}".format( xmax ),
          "xtitle         {}".format( xtitle ),
          "ymin           0",
          "subymin        {}".format( subymin ),
          "subymax        2",
          "rebin          {}".format( rebin ),
          "plotImpact     false"
        ]

  lines = "\n".join(lines)+"\n"

  file = open( "plot_pars.txt", "w+" )
  file.write( lines )
  file.close()

if __name__ == "__main__":
  main()
