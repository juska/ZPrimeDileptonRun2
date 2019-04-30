# ZDileptonAnalysis2017

1==> Setting up analyzer workflow and libraries for MT2
=======================================================

mkdir 2017DATA_Zprime

cd 2017DATA_Zprime

cmsrel CMSSW_9_4_12

cd CMSSW_9_4_12/src

cmsenv

cd ..

mkdir oxbridge

cd /uscms_data/d3/username

mkdir oxbridge

cd oxbridge

wget http://www.hep.phy.cam.ac.uk/~lester/dtm662/mt2/Releases/oxbridgekinetics.tar.gz

tar -xzvf oxbridgekinetics.tar.gz

cd oxbridgekinetics-1.2

./configure --prefix=/uscms_data/d3/username/2017DATA_Zprime/CMSSW_9_4_12/oxbridge/install

(still inside oxbridgekinetics-1.2)

make

make install

then in

/uscms_data/d3/username/2017DATA_Zprime/CMSSW_9_4_12/config/toolbox/slc6_amd64_gcc630/tools/selected

copy:

cp /uscms_data/d3/cihar29/newAnalysis/CMSSW_9_4_12/config/toolbox/slc6_amd64_gcc630/tools/selected/oxbridgekinetics-1.2.xml .

then in CMSSW_9_4_12

scram setup oxbridgekinetics-1.2

scram b

then in src

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

2==> analyzer workflow
=======================================================
