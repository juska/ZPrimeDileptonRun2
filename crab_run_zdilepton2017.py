import time
from WMCore.Configuration import Configuration
config = Configuration()

isMC = True

if isMC:
  inputFiles = []
  outputFiles = ["analysis_MC.root"]
  splitting = 'FileBased'
  lumiMask = ''
  unitsPerJob = 10

else:
  inputFiles = ["pileup_json.txt","json.txt"]
  outputFiles = ["analysis_Data.root"]
  splitting = 'LumiBased'
  lumiMask = 'json.txt'
  unitsPerJob = 100

config.section_("General")
config.General.requestName = ''
config.General.workArea = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/ConfFile_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = inputFiles
config.JobType.outputFiles = outputFiles

config.section_("Data")
config.Data.inputDataset = ''

 #'/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD'
 #'/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD'
 #'/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD'
 #'/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD'
 #'/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD'

config.Data.splitting = splitting
config.Data.lumiMask = lumiMask

config.Data.unitsPerJob = unitsPerJob
config.Data.outLFNDirBase = '/store/user/broozbah/'
config.Data.publication = False
#config.Data.ignoreLocality = True
#config.Data.publishDataName = 'analysis_tree'

config.section_("Site")
#config.Site.blacklist = ['T1_US_FNAL']
#config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = "T3_US_FNALLPC"


if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    sample = {}
    version = {}

    f = open( "RSGluonToTT.txt", "r" )
    datasets = []
    for line in f:
        datasets.append(line.rstrip())

    for dataset in datasets:
        sample.setdefault(dataset, [])
        sample[dataset].append(dataset.split("/")[1])
        version.setdefault(dataset, [])
        version[dataset].append(dataset.split("/")[2].split("_realistic")[1])
    for key in version:
        config.General.requestName = sample[key][0]+version[key][0]
        config.Data.inputDataset = key
        config.Data.outLFNDirBase = '/store/user/broozbah/ZPRIME_2017/Signal/' + sample[key][0]+version[key][0]
        crabCommand('submit', config = config) 
        time.sleep(40)   

# source /cvmfs/cms.cern.ch/crab3/crab.csh
