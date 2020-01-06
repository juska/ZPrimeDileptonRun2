#chad harrington 2019
#example: python python/run_all.py -drs OFF ON -c mm -uncerts jec pileup muIdSys --makePileupRewgtFile
#run everything: python -u python/run_all.py > & log.txt &
#simple tests: python -u python/run_all.py -drs ON -c mm --noUncerts > & log.txt &
import os
import argparse
import time
from files import mcTups, dataDict
from calculate_pileup import calc_pileup
from calculate_weights import calc

drs_default = ["OFF", "ON", "ONbt", "ONnb", "REVERSE"]
channels_default = ["mm", "em", "ee"]
uncerts_default = ["jec", "jer", "btagSF", "mistagSF", "pileup", "topPtWeight", "pdf", "q2ttbar", "muTrigSys", "muIdSys", "eleTrigSys", "eleIdSys"]
types = ["UP", "DOWN"]

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
  parser.add_argument( "-uncerts", "--uncertainties", nargs='+', required=False,
                       help="uncertainties to analyze, like -uncerts jec pileup muIdSys",
                       default=uncerts_default
                     )
  parser.add_argument( "--makePileupRewgtFile", action='store_true', required=False, help="If enabled, create the pileup reweighting root file" )
  parser.add_argument( "--makeBtagFile",        action='store_true', required=False, help="If enabled, create the btag efficiency root file" )
  parser.add_argument( "--noUncerts",           action='store_true', required=False, help="If enabled, do not run over any uncertainties" )
  return parser.parse_args()

def main() :
  start = time.time()

  args = handleArgs()
  if args.makePileupRewgtFile :
    calc_pileup()
  if args.makeBtagFile :
    os.system( "root -l -b -q btag_efficiency.c" )

  drs = args.drCuts
  channels = args.channels
  uncerts = args.uncertainties
  if args.noUncerts : uncerts = []

  print "drCuts:", " ".join(drs)
  print "channels:", " ".join(channels)
  print "uncertainties:", " ".join(uncerts)

  era = "Fall17_17Nov2017"
  periods = ("B", "C", "DE", "F")
  version = "32"
  eras = "" 
  for p in periods :
    eras += era + p + "_V" + version + "_DATA "

  for drcut in drs :
    if drcut not in drs_default : print drcut + " not found!"; continue
    if os.path.isdir( drcut ) : raise Exception(drcut + " directory already exists!")
    os.system( "mkdir " + drcut )
    os.system( "mkdir " + drcut + "/logs" )

    for channel in channels :
      if channel not in channels_default : print channel + " not found!"; continue
      os.system( "mkdir " + drcut + "/" + channel )

      dataTup = dataDict[channel]
      print "Processing:", drcut, channel, dataTup.name
      writePars( "false", dataTup.dir, dataTup.name, drcut, channel, eras )
      file = open( drcut + "/logs/log_" + channel + ".txt", "w+" )
      file.write( os.popen( "analyzer pars.txt" ).read() )
      file.flush()

      calc( dataTup.lumi )

      for mcTup in mcTups :
        print "Processing:", drcut, channel, mcTup.name
        writePars( "true", mcTup.dir, mcTup.name, drcut, channel, era+"_V"+version+"_MC" )
        file.write( os.popen( "analyzer pars.txt mc_weights.txt" ).read() )
        file.flush()
      file.close()

      for uncert in uncerts :
        if uncert not in uncerts_default : print uncert + " not found!"; continue
        for type in types :

          file = open( drcut + "/logs/log_" + channel + "_" + uncert+type + ".txt", "w+" )
          for mcTup in mcTups :
            print "Processing:", drcut, channel, uncert+type, mcTup.name
            writePars( "true", mcTup.dir, mcTup.name, drcut, channel, era+"_V"+version+"_MC", uncert+type )
            file.write( os.popen( "analyzer pars.txt mc_weights.txt" ).read() )
            file.flush()
          file.close()

  print "{} seconds".format( time.time()-start )

def writePars( isMC, inDir, fname, drcut, channel, eras, uncert="" ) :

  lines=[ "isMC           %s" % isMC,
          "inName         %s/%s.root" % (inDir, fname),
          "outName        %s/%s/%s_%s.root" % (drcut, channel, fname, channel),
          "drCut          %s" % drcut,
          "channel        %s" % channel,
          "eras           %s" % eras,
          "jet_type       AK4PFchs"
        ]

  if isMC == "true" :
    lines.append( "res_era        Fall17_V3_MC" )
    lines.append( "muTrigSfName   /uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root" )
    lines.append( "muIdSfName     /uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/RunBCDEF_SF_ID_syst.root" )
    lines.append( "eTrigSfName    /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/electronTrigSF.root" )
    lines.append( "eRecoSfName    /uscms_data/d3/broozbah/ZPRIME_2017/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root" )
    lines.append( "eIdSfName      /uscms_data/d3/broozbah/ZPRIME_2017/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/2017_ElectronTight.root" )
    lines.append( "btagName       /uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/btag_eff.root" )
    lines.append( "pileupName     /uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/mu_weights.root" )
  if uncert != "" :
    lines.append( "uncert         %s" % uncert )

  lines = "\n".join(lines)+"\n"

  file = open( "pars.txt", "w+" )
  file.write( lines )
  file.close()

if __name__ == "__main__":
    main()
