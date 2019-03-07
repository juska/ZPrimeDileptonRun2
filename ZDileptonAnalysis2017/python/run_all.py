import os

def main() :

  drs = ("OFF", "ON", "ONbt", "ONnb", "REVERSE")
  channels = ("mm", "em", "ee")

  uncerts = ("jec", "jer", "btagSF", "mistagSF", "pileup", "topPtWeight", "pdf", "q2ttbar", "muTrigSys", "muIdSys", "eleTrigSys", "eleIdSys")
  types = ("UP", "DOWN")

  cdir = "root://131.225.204.161:1094//store/user/cihar29/2017Analysis_trees"
  bdir = ""
  mcSets = ( (cdir,"TTjets"), (cdir,"dylow") )

  era = "Fall17_17Nov2017"
  periods = ("B", "C", "DE", "F")
  version = "32"
  eras = "" 
  for p in periods :
    eras += era + p + "_V" + version + "_DATA "

  #run pileup_weights (pass datasets from here (make exe in bin))
  #make timer

  for dr in drs :
#    os.system( "mkdir " + dr )
    for channel in channels :
#      os.system( "mkdir " + dr + "/" + channel )
#      os.system( "mkdir " + dr "/logs" )

      dir, fname, lumi = cdir, "Muon", "41.86"
      if channel == "ee" :
        dir, fname, lumi = bdir, "Ele", "41.53"

      writePars( "false", dir, fname, dr, channel, eras )
#      os.system( "analyzer pars.txt > logs/log_" + set + "_" + channel + ".txt" )

      os.system( "python python/calculate_weights.py -l " + lumi )

      for set in mcSets :
        writePars( "true", set[0], set[1], dr, channel, era+"_V"+version+"_MC" )
#        os.system( "analyzer pars.txt mc_weights.txt > logs/log_" + set + "_" + channel + ".txt" )

        for uncert in uncerts :
          for type in types :
            writePars( "true", set[0], set[1], dr, channel, era+"_V"+version+"_MC", uncert+type )
#            os.system( "analyzer pars.txt mc_weights.txt > logs/log_" + set + "_" + channel + "_" + uncert+type + ".txt" )

#   cat all logs together

def writePars( isMC, inDir, fname, dr, channel, eras, uncert="" ) :

  lines=[ "isMC           %s" % isMC,
          "inName         %s/%s.root" % (inDir, fname),
          "outName        %s/%s/%s_%s.root" % (dr, channel, fname, channel),
          "drCut          %s" % dr,
          "channel        %s" % channel,
          "eras           %s" % eras,
          "jet_type       AK4PFchs"
        ]

  if isMC == "true" :
    lines.append( "res_era        Fall17_V3_MC" )
    lines.append( "muTrigSfName   /uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root" )
    lines.append( "muIdSfName     /uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/src/ZDileptonAnalysis2017/ZDileptonAnalysis2017/RunBCDEF_SF_ID_syst.root" )
    lines.append( "eTrigSfName    eTrigSF" )
    lines.append( "eRecoSfName    eRecoSF" )
    lines.append( "eIdSfName      eIdSF" )
    lines.append( "btagName       btag" )
    lines.append( "pileupName     pileup" )
  if uncert != "" :
    lines.append( "uncert         %s" % uncert )

  lines = "\n".join(lines)

  file = open( "pars.txt", "w+" )
  file.write( lines )
  file.close()

if __name__ == "__main__":
    main()
