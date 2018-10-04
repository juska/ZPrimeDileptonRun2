void counter(){

  TFile* dataFile = TFile::Open("analysis_Data.root");
  TFile* mcFile = TFile::Open("analysis_MC.root");

  vector<string> sample = {"Data","MC"};
  enum Cuts{
    countEvts, countDilep, countLeppt, countDilepmass, countJetpteta, countMet, numCuts
  };
  map <string, vector<pair<string, double> > > v_cuts;
  v_cuts["Data"]= { make_pair("Initial",0.),make_pair("Dilepton selection",0.),
                    make_pair("Lepton Pt Cut",0.),make_pair("Dilepton Mass Cut",0.),
                    make_pair("Leading Jet Pt/eta cut",0.),make_pair("MET Filters",0.)};
  v_cuts["MC"]= { make_pair("Initial",0.),make_pair("Dilepton selection",0.),
                    make_pair("Lepton Pt Cut",0.),make_pair("Dilepton Mass Cut",0.),
                    make_pair("Leading Jet Pt/eta cut",0.),make_pair("MET Filters",0.)};
  double weight0 = 1;
  TIter datakey(dataFile->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey*)datakey()) ) {
    TString keyname = key->GetName();
    if      (keyname=="totalEvts")     v_cuts["Data"][countEvts].second      += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="dilep_cut")     v_cuts["Data"][countDilep].second     += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="leppt_cut")     v_cuts["Data"][countLeppt].second     += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="dilepmass_cut") v_cuts["Data"][countDilepmass].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="jetpteta_cut")  v_cuts["Data"][countJetpteta].second  += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="met_cut")       v_cuts["Data"][countMet].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
  }
  TIter mckey(mcFile->GetListOfKeys());
  while ( (key = (TKey*)mckey()) ) {
    TString keyname = key->GetName();
    if      (keyname=="totalEvts")     v_cuts["MC"][countEvts].second      += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="dilep_cut")     v_cuts["MC"][countDilep].second     += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="leppt_cut")     v_cuts["MC"][countLeppt].second     += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="dilepmass_cut") v_cuts["MC"][countDilepmass].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="jetpteta_cut")  v_cuts["MC"][countJetpteta].second  += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="met_cut")       v_cuts["MC"][countMet].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
  }
  printf("Data file: /store/data/Run2017B/SingleMuon/MINIAOD/17Nov2017-v1/70000/FEDE478C-D3D7-E711-A54E-02163E011BE4.root\n");
  printf("===============================================================================================================\n");
  printf("||%5s%-30s%5s||%5s%-10s%5s||%10s%-20s%3s||\n", "","","", "","Nevents","", "","Efficiency (Relative efficiency)","");
  for (int i=0; i<numCuts; i++){
    printf("||%5s%-30s%5s||%5s%-10.1f%5s||%10s%-14.6f (%-4.6f)%10s||\n", "",v_cuts["Data"][i].first.c_str(),"", "",v_cuts["Data"][i].second,"", "",i==countEvts ? 1 : v_cuts["Data"][i].second/v_cuts["Data"][0].second, i==countEvts ? 1 : v_cuts["Data"][i].second/v_cuts["Data"][i-1].second,"");
  }


  printf("MC file: /store/mc/RunIIFall17MiniAODv2/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/00A57BF8-DD41-E811-82C4-008CFAED6D70.root\n");
  printf("===============================================================================================================\n");
  printf("||%5s%-30s%5s||%5s%-10s%5s||%10s%-20s%3s||\n", "","","", "","Nevents","", "","Efficiency (Relative efficiency)","");
  for (int i=0; i<numCuts; i++){
    printf("||%5s%-30s%5s||%5s%-10.1f%5s||%10s%-14.6f (%-4.6f)%10s||\n", "",v_cuts["MC"][i].first.c_str(),"", "",v_cuts["MC"][i].second,"", "",i==countEvts ? 1 : v_cuts["MC"][i].second/v_cuts["MC"][0].second, i==countEvts ? 1 : v_cuts["MC"][i].second/v_cuts["MC"][i-1].second,"");
  }

}
