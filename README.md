# ZPrimeDileptonRun2

1==> Setting up analyzer workflow
=======================================================


# Best CMSSW version to be checked

cmsrel CMSSW_10_6_8

cd CMSSW_10_6_8/src && cmsenv

scram b -j 8

git clone https://github.com/juska/ZPrimeDileptonRun2

(copy all related libraries for analysis (can copy from /uscms_data/d3/broozbah/ZPRIME_2017/CMSSW_9_4_12/src/))

=======================================================

2==> Skimmer workflow
=======================================================

source file: ZDileptonAnalysis2017/ZDileptonAnalysis2017/plugins/ZDileptonAnalysis2017.cc

Configuration file: DileptonAnalysis2017/ZDileptonAnalysis2017/python/ConfFile_cfg.py
(flag to be modified is isMC)

Crab configuration file to submit multiple samples at once: ZDileptonAnalysis2017/ZDileptonAnalysis2017/crab_run_zdilepton2017.py

to run it, input your sample/samples in a text file and input the text file in crab_run_zdilepton2017.py, line 62 (e.g. RSGluonToTT.txt)

command for running is "python crab_run_zdilepton2017.py" 

=======================================================

3==> analyzer workflow
=======================================================

Analyzer source file : ZDileptonAnalysis2017/ZDileptonAnalysis2017/bin/analyzer.cc

Analyzer runner :  ZDileptonAnalysis2017/ZDileptonAnalysis2017/python/run_all.py 

Instruction on the selections while running analyzer is documented:

https://github.com/broozbah/ZDileptonAnalysis2017/blob/master/ZDileptonAnalysis2017/python/run_all.py#L2:L4


=======================================================

4==> plotter workflow
=======================================================

Analyzer source file : ZDileptonAnalysis2017/ZDileptonAnalysis2017/bin/plotter.cc

Analyzer runner :  ZDileptonAnalysis2017/ZDileptonAnalysis2017/python/plot_all.py 

Instruction on the selections while running analyzer is documented:

https://github.com/broozbah/ZDileptonAnalysis2017/blob/master/ZDileptonAnalysis2017/python/plot_all.py#L2


=======================================================

5==> createLog workflow
=======================================================


Analyzer source file : ZDileptonAnalysis2017/ZDileptonAnalysis2017/bin/createLog.cc

Analyzer runner :  executable createLog 

Instruction on the selections while running analyzer is documented:

https://github.com/broozbah/ZDileptonAnalysis2017/blob/master/ZDileptonAnalysis2017/bin/createLog.cc#L1:L3
