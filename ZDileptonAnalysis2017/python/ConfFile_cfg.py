import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [
   #'/store/data/Run2017B/SingleMuon/MINIAOD/17Nov2017-v1/70000/FEDE478C-D3D7-E711-A54E-02163E011BE4.root'
   '/store/data/Run2017B/DoubleEG/MINIAOD/17Nov2017-v1/70000/A8B1289A-91E0-E711-BBD2-F4E9D4AF7940.root'
] );

isMC = cms.bool(False)
gts = {'BCDEFG':'80X_dataRun2_2016SeptRepro_v7', 'H':'80X_dataRun2_Prompt_v16'}

if isMC:
  OutputName = "_MC"  
  metLabel = "SIM"
  gt = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

else:
  OutputName = "_Data"
  metLabel = "RECO"
  gt = gts['H']


dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

my_eid_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
for idmod in my_eid_modules:
setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.demo = cms.EDAnalyzer('ZDileptonAnalysis2017',
    fileName = cms.string("analysis" + OutputName + ".root"),
    isMC = isMC,
    genParticleTag = cms.InputTag("prunedGenParticles"),
    muTag = cms.InputTag("slimmedAddPileupInfo"),
    muonTag = cms.InputTag("slimmedMuons"),
    electronTag = cms.InputTag("slimmedElectrons"),
    eleVetoIdMapToken = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto")
)


process.myseq = cms.Sequence( process.demo )

process.p = cms.Path(process.myseq)
