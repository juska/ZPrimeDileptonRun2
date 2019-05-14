# ZDileptonAnalysis2017

1==> Setting up analyzer workflow and libraries for MT2
=======================================================

cmsrel CMSSW_9_4_12

cd CMSSW_9_4_12/ && cmsenv && mkdir oxbridge && cd ..

wget http://www.hep.phy.cam.ac.uk/~lester/dtm662/mt2/Releases/oxbridgekinetics.tar.gz

tar -xzvf oxbridgekinetics.tar.gz

cd oxbridgekinetics-1.3

./configure --prefix=$PWD/../CMSSW_9_4_12/oxbridge/install

make

make install

cd ../CMSSW_9_4_12/

wget https://raw.githubusercontent.com/broozbah/ZDileptonAnalysis2017/master/oxbridgekinetics-1.3.xml -P config/toolbox/slc6_amd64_gcc630/tools/selected/

scram setup oxbridgekinetics-1.3

scram b -j 8

cd src

git clone https://github.com/broozbah/ZDileptonAnalysis2017

copy all related libraries for analysis (can copy from /uscms_data/d3/broozbah/ZPRIME_2017/CMSSW_9_4_12/src/)

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
