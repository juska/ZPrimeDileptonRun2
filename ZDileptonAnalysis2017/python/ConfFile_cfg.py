import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [
   #'/store/data/Run2017B/SingleMuon/MINIAOD/17Nov2017-v1/70000/FEDE478C-D3D7-E711-A54E-02163E011BE4.root'
   #'/store/data/Run2017B/DoubleEG/MINIAOD/17Nov2017-v1/70000/A8B1289A-91E0-E711-BBD2-F4E9D4AF7940.root'
   '/store/mc/RunIIFall17MiniAODv2/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/00A57BF8-DD41-E811-82C4-008CFAED6D70.root'
] );

isMC = cms.bool(True)
#gts = {'BCDEFG':'80X_dataRun2_2016SeptRepro_v7', 'H':'80X_dataRun2_Prompt_v16'}

if isMC:
  OutputName = "_MC"  
  metLabel = "SIM"
  gt = '94X_mc2017_realistic_v14'
  inData = False

else:
  OutputName = "_Data"
  metLabel = "RECO"
  gt = '94X_dataRun2_ReReco17_forValidation'
  inData = True

process.load( "Configuration.Geometry.GeometryIdeal_cff" )
process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, gt)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD (
    process,
    isData = inData,
    fixEE2017 = True,
    fixEE2017Params = {'userawPt': True, 'PtThreshold':50.0, 'MinEtaThreshold':2.65, 'MaxEtaThreshold': 3.139} ,
    postfix = "ModifiedMET"
)

#process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
#process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
#process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
#process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
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
    triggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    prescalesTag = cms.InputTag("patTrigger")   
)


process.myseq = cms.Sequence( process.fullPatMetSequenceModifiedMET * process.egmGsfElectronIDSequence * process.demo )

process.p = cms.Path(process.myseq)
