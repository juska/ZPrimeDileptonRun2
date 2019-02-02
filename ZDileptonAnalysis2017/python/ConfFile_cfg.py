import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [
   #'/store/data/Run2017B/SingleMuon/MINIAOD/17Nov2017-v1/70000/FEDE478C-D3D7-E711-A54E-02163E011BE4.root'
   #'/store/data/Run2017B/DoubleEG/MINIAOD/17Nov2017-v1/70000/A8B1289A-91E0-E711-BBD2-F4E9D4AF7940.root'
   #'/store/mc/RunIIFall17MiniAODv2/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/00A57BF8-DD41-E811-82C4-008CFAED6D70.root'
   '/store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/100000/001642F1-6638-E811-B4FA-0025905B857A.root'
   #'/store/data/Run2017B/DoubleEG/MINIAOD/31Mar2018-v1/30000/0EF2F0D3-9137-E811-B356-6CC2173DA930.root'
] );

isMC = cms.bool(False)

if isMC:
  OutputName = "_MC"  
  metLabel = "SIM"
  gt = '94X_mc2017_realistic_v14'

else:
  OutputName = "_Data"
  metLabel = "RECO"
  gt = '94X_dataRun2_ReReco_EOY17_v6'

process.load( "Configuration.Geometry.GeometryIdeal_cff" )
process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, gt)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD (
    process,
    isData = not isMC,
    fixEE2017 = True,
    fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
    postfix = "ModifiedMET"
)

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )

#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#setupEgammaPostRecoSeq(process,
#                       era='2017-Nov17ReReco',
#                       runVID=True)

switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_eid_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']
for idmod in my_eid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.demo = cms.EDAnalyzer('ZDileptonAnalysis2017',
    fileName = cms.string("analysis" + OutputName + ".root"),
    isMC = isMC,
    genParticleTag = cms.InputTag("prunedGenParticles"),
    genEventTag = cms.InputTag("generator"),
    extLHETag = cms.InputTag("externalLHEProducer"),
    muTag = cms.InputTag("slimmedAddPileupInfo"),
    muonTag = cms.InputTag("slimmedMuons"),
    electronTag = cms.InputTag("slimmedElectrons"),
    minLepPt = cms.double(45.),
    minSubLepPt = cms.double(25.),
    minDiLepMass = cms.double(20.),
    jetTag = cms.InputTag("slimmedJets"),
    puppiJetTag = cms.InputTag("slimmedJetsPuppi"),
    btag = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    btag_deepcsvprobb = cms.string("pfDeepCSVJetTags:probb"),
    btag_deepcsvprobbb = cms.string("pfDeepCSVJetTags:probbb"),
    minLeadJetPt = cms.double(80.),
    patTrgLabel = cms.InputTag("TriggerResults", "", "RECO"),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    genJetTag = cms.InputTag("slimmedGenJets"),
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),
    effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
    metTag = cms.InputTag("slimmedMETs"),
    puppiMetTag = cms.InputTag("slimmedMETsPuppi"),
    triggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    prescalesTag = cms.InputTag("patTrigger")   
)

process.myseq = cms.Sequence( process.fullPatMetSequenceModifiedMET *
                              #process.egammaPostRecoSeq *
                              process.egmGsfElectronIDSequence *
                              process.ecalBadCalibReducedMINIAODFilter *
                              process.demo )

process.p = cms.Path(process.myseq)
