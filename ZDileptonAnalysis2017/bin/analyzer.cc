//charles harrington and bahareh roozbahani
//execute as analyze pars.txt mc_weights.txt

#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"

#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <time.h>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/Math/interface/deltaPhi.h"

using namespace std;

void FillHist1D(const TString& histName, const Double_t& value, const double& weight);
void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight);
void setPars(const string& parFile);
void setWeight(const string& parFile);
bool sortJetPt(const pair<int, float>& jet1, const pair<int, float>& jet2){ return jet1.second > jet2.second; }
bool newBTag( const float& coin, const float& pT, const int& flavor, const bool& oldBTag, TGraphAsymmErrors& g_eff, const TString& variation );
double rms_pm(const vector<float>& vec);

map<TString, TH1*> m_Histos1D;
map<TString, TH2*> m_Histos2D;

//parameters- edit in pars.txt
bool isMC;
TString jec="NOMINAL", jer="NOMINAL", btagSF="NOMINAL", mistagSF="NOMINAL", pileup="NOMINAL", pdf="NOMINAL", q2ttbar="NOMINAL"; //NOMINAL, UP, DOWN
TString muTrigSys="NOMINAL", muIdSys="NOMINAL", eleTrigSys="NOMINAL", eleIdSys="NOMINAL";
TString topPtWeight="NOMINAL"; //NOMINAL (sqrt tPt*tbarPt), UP (no top reweighting), DOWN is NOT an option (made symmetric with UP in later stages)
TString setSUMDRCut="OFF";
//ON (keep events with sumrmin<2.0), REVERSE (keep events with sumrmin>=2.0), OFF (no cut)
//ONbt (keep events with sumrmin<1.0),  ONnb (keep events with 1<=sumrmin<2.0)
TString inName, outName, muTrigSfName, muTrackSfName, eTrigSfName, eRecoSfName, eIdSfName, btagName, pileupName;
map<TString, double> muTrigSfNames, muIdSfNames;
string channel, jet_type, res_era;
vector<string> eras;
double weight0, weight;

const int MAXJET = 50;
const int MAXLEP = 20;
const int MAXGEN = 20;
const float MUONMASS = 0.10566;
const float ELEMASS = 0.;
const float btagWP_L = 0.5426;
const float btagWP_M = 0.8484;

int main(int argc, char* argv[]){

  if (argc == 1) { cout << "Please provide a parameter file." << endl; return -1; }

  string parFile = argv[1];
  setPars(parFile);

  weight0 = -1;
  if (argc == 3)     { string wFile = argv[2]; setWeight(wFile); }
  if (weight0 == -1) { cout << "Weight set to 1" << endl; weight0 = 1.; }
  else                 cout << "Weight set to " << weight0 << endl;

  //Jet Corrections and Resolution//

  /*map<string, JetCorrectorParameters*> ResJetPars, L3JetPars, L2JetPars, L1JetPars;
  map<string, vector<JetCorrectorParameters> > jetPars, jetL1Pars;
  map<string, FactorizedJetCorrector*> jetCorrectors, jetL1Correctors;
  map<string, JetCorrectionUncertainty*> jecUncert;
  map<string, pair<int, int> > m_IOV;

  cout << endl << "Using eras: " << endl;
  for(vector<string>::iterator i_era = eras.begin(); i_era != eras.end(); ++i_era) {
    string era = *i_era;

    ResJetPars[era] = new JetCorrectorParameters(era + "/" + era + "_L2L3Residual_" + jet_type + ".txt");
    L3JetPars[era] = new JetCorrectorParameters(era + "/" + era + "_L3Absolute_" + jet_type + ".txt");
    L2JetPars[era] = new JetCorrectorParameters(era + "/" + era + "_L2Relative_" + jet_type + ".txt");
    L1JetPars[era] = new JetCorrectorParameters(era + "/" + era + "_L1FastJet_" + jet_type + ".txt");
    jecUncert[era] = new JetCorrectionUncertainty(era + "/" + era + "_Uncertainty_" + jet_type + ".txt");

    jetPars[era].push_back( *L1JetPars[era] );
    jetPars[era].push_back( *L2JetPars[era] );
    jetPars[era].push_back( *L3JetPars[era] );
    jetPars[era].push_back( *ResJetPars[era] );

    jetCorrectors[era] = new FactorizedJetCorrector( jetPars[era] );

    jetL1Pars[era].push_back( *L1JetPars[era] );
    jetL1Correctors[era] = new FactorizedJetCorrector( jetL1Pars[era] );

    if (isMC) cout << era << endl;
    else {
      TString tera = era.data();

      if ( tera.Contains("BCDV", TString::kIgnoreCase) ) { cout << "BCD "; m_IOV[era] = make_pair(1, 276811); }
      else if ( tera.Contains("EFV", TString::kIgnoreCase) ) { cout << "EF "; m_IOV[era] = make_pair(276831, 278801); }
      else if ( tera.Contains("GV", TString::kIgnoreCase) ) { cout << "G "; m_IOV[era] = make_pair(278802, 280385); }
      else if ( tera.Contains("HV", TString::kIgnoreCase) ) { cout << "H "; m_IOV[era] = make_pair(280919, 300000); }

      cout << "[" << m_IOV[era].first << ", " << m_IOV[era].second << "]" << endl;
    }
  }

  JME::JetResolution res_obj;
  JME::JetResolutionScaleFactor ressf_obj;

  if (isMC) {
    res_obj = JME::JetResolution( res_era + "/" + res_era + "_PtResolution_" + jet_type + ".txt" );
    ressf_obj = JME::JetResolutionScaleFactor( res_era + "/" + res_era + "_SF_" + jet_type + ".txt" );
  }
  */
  //Open Files//

  TFile* inFile = TFile::Open(inName);

  TTree* T = (TTree*) inFile->Get("T");
  Long64_t nEntries = T->GetEntries();
  cout << endl << nEntries << " Events" << endl;
  cout << "Processing " + inName << endl;
  TString name = inName( inName.Last('/')+1, inName.Last('.')-inName.Last('/')-1);
  //cout << "Channel: " + channel << endl;

  //Reweighting and SF Files//
  /*
  TH1F* pileup_weights=0;
  TH2F* muTrigSfHist=0, *muIdSfHist=0, *eTrigSfHist=0, *eRecoSfHist=0, *eIdSfHist=0;
  TGraphAsymmErrors* muTrackSfGraph=0, *btag_eff_b=0, *btag_eff_c=0, *btag_eff_udsg=0;

  float muTrig_pT=0, muId_pT=0, eTrig_pT=0, eReco_pT=0, eId_pT=0;

  if (isMC) {
    vector<TH2F*> muTrigSfHists;
    for (const auto& it : muTrigSfNames) {
      TFile* muTrigSfFile = TFile::Open(it.first);
      muTrigSfHists.push_back( (TH2F*) muTrigSfFile->Get("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio") );
      muTrigSfHists.back()->Scale(it.second);
    }
    muTrigSfHist = muTrigSfHists[0];
    for (unsigned int i=1; i<muTrigSfHists.size(); i++) muTrigSfHist->Add( muTrigSfHists[i] );

    vector<TH2F*> muIdSfHists;
    for (const auto& it : muIdSfNames) {
      TFile* muIdSfFile = TFile::Open(it.first);
      muIdSfHists.push_back( (TH2F*) muIdSfFile->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio") );
      muIdSfHists.back()->Scale(it.second);
    }
    muIdSfHist = muIdSfHists[0];
    for (unsigned int i=1; i<muIdSfHists.size(); i++) muIdSfHist->Add( muIdSfHists[i] );

    TFile* muTrackSfFile = TFile::Open(muTrackSfName);
    muTrackSfGraph = (TGraphAsymmErrors*) muTrackSfFile->Get("ratio_eff_eta3_dr030e030_corr");

    TFile* eTrigSfFile = TFile::Open(eTrigSfName);
    eTrigSfHist = (TH2F*) eTrigSfFile->Get("tight_ScaleFactor");

    TFile* eRecoSfFile = TFile::Open(eRecoSfName);
    eRecoSfHist = (TH2F*) eRecoSfFile->Get("EGamma_SF2D");

    TFile* eIdSfFile = TFile::Open(eIdSfName);
    eIdSfHist = (TH2F*) eIdSfFile->Get("EGamma_SF2D");

    muTrig_pT = muTrigSfHist->GetYaxis()->GetBinCenter(muTrigSfHist->GetYaxis()->GetNbins());
    muId_pT = muIdSfHist->GetYaxis()->GetBinCenter(muIdSfHist->GetYaxis()->GetNbins());
    eTrig_pT = eTrigSfHist->GetYaxis()->GetBinLowEdge(eTrigSfHist->GetYaxis()->GetNbins());
    eReco_pT = eRecoSfHist->GetYaxis()->GetBinCenter(eRecoSfHist->GetYaxis()->GetNbins());
    eId_pT = eIdSfHist->GetYaxis()->GetBinCenter(eIdSfHist->GetYaxis()->GetNbins());

    TFile* btagFile = TFile::Open(btagName);
    btag_eff_b = (TGraphAsymmErrors*) btagFile->Get( Form("%s/b_LWP_%s", channel.data(), channel.data()) );
    btag_eff_c = (TGraphAsymmErrors*) btagFile->Get( Form("%s/c_LWP_%s", channel.data(), channel.data()) );
    btag_eff_udsg = (TGraphAsymmErrors*) btagFile->Get( Form("%s/udsg_LWP_%s", channel.data(), channel.data()) );

    TFile* pileupFile = TFile::Open(pileupName);

    TIter nextHist(pileupFile->GetDirectory(pileup)->GetListOfKeys());
    TKey* histKey;
    TString name = inName( inName.Last('/')+1, inName.Index('.')-inName.Last('/')-1 );
    while ( (histKey = (TKey*)nextHist()) ) {
      TString keyname = histKey->GetName();

      if ( name.EqualTo( keyname(0, keyname.Last('_')), TString::kIgnoreCase ) ) {
        pileup_weights = (TH1F*) histKey->ReadObj();
        break;
      }
    }
  }

  //Skims and Cuts//

  enum Cuts{
    countEvts, countDilep, countLeppt, countDilepmass, countJetpteta, countMet,
    channelCut, trigCut, lepkinCut, signCut, thirdLepCut, dilepmassCut, dilepVetoCut, ptrelCut, dRCut, jetCut,
    zerobtagCut1jet_metR, onebtagCut1jet_metR, zerobtagCut2jets_metR, morethan0btagCut2jets_metR,
    zerobtagCut1jet, onebtagCut1jet, zerobtagCut2jets, morethan0btagCut2jets, numCuts
  };
  vector<pair<string, double> > v_cuts(numCuts);

  v_cuts[countEvts]=make_pair("Initial",0.); v_cuts[countDilep]=make_pair("Dilepton selection",0.); v_cuts[countLeppt]=make_pair("Lepton Pt Cut",0.);
  v_cuts[countDilepmass]=make_pair("Dilepton Mass Cut",0.); v_cuts[countJetpteta]=make_pair("Leading Jet Pt/eta cut",0.);
  v_cuts[countMet]=make_pair("MET Filters",0.); v_cuts[channelCut]=make_pair("Correct Channel",0.); v_cuts[signCut]=make_pair("Opposite Lepton Sign",0.);
  v_cuts[trigCut]=make_pair("HLT Trigger",0.); v_cuts[lepkinCut]=make_pair("Lepton kinematics cut",0.); v_cuts[thirdLepCut]=make_pair("Third lepton cut",0.);
  v_cuts[dilepmassCut]=make_pair("Dilepton mass cut",0.); v_cuts[ptrelCut]=make_pair("pTrel cut",0.); v_cuts[dRCut]=make_pair("DeltaR cut",0.);
  v_cuts[dilepVetoCut]=make_pair("Z-mass veto",0.); v_cuts[jetCut]=make_pair(">= 1 jet",0.);
  v_cuts[zerobtagCut1jet_metR]=make_pair("= 1 Jet, = 0 btags, metR",0.);  v_cuts[onebtagCut1jet_metR]=make_pair("= 1 Jet, >= 1 btag, metR",0.);
  v_cuts[zerobtagCut2jets_metR]=make_pair(">= 2 Jets, = 0 btags, metR",0.); v_cuts[morethan0btagCut2jets_metR]=make_pair(">= 2 Jets, >= 1 btag, metR",0.);
  v_cuts[zerobtagCut1jet]=make_pair("= 1 Jet, = 0 btags",0.);  v_cuts[onebtagCut1jet]=make_pair("= 1 Jet, >= 1 btag",0.);
  v_cuts[zerobtagCut2jets]=make_pair(">= 2 Jets, = 0 btags",0.); v_cuts[morethan0btagCut2jets]=make_pair(">= 2 Jets, >= 1 btag",0.);

  TIter nextkey(inFile->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey*)nextkey()) ) {
    TString keyname = key->GetName();

    if      (keyname=="totalEvts")     v_cuts[countEvts].second      += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="dilep_cut")     v_cuts[countDilep].second     += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="leppt_cut")     v_cuts[countLeppt].second     += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="dilepmass_cut") v_cuts[countDilepmass].second += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="jetpteta_cut")  v_cuts[countJetpteta].second  += weight0 * (*(vector<int>*)key->ReadObj())[0];
    else if (keyname=="met_cut")       v_cuts[countMet].second       += weight0 * (*(vector<int>*)key->ReadObj())[0];
    //else if (keyname=="filter_failed")
  }
  if ( (int) (v_cuts[countMet].second + 0.5) != (int) (weight0 * nEntries + 0.5) ) { cout << "hadd added incorrectly." << endl; return -1; }

  vector<pair<string, double> > v_llcuts = v_cuts, v_ljcuts = v_cuts, v_jjcuts = v_cuts;
  vector<pair<string, double> >* v_cuts_ptr = &v_llcuts;

  //ttbar reweighting
  if ( inName.Contains("ttbar", TString::kIgnoreCase) && topPtWeight=="NOMINAL" ) {
    double topPtWeightNOM=0, total=0;

    nextkey = inFile->GetListOfKeys();
    while ( (key = (TKey*)nextkey()) ) {
      TString keyname = key->GetName();

      if      (keyname=="totalEvts")      total          += (*(vector<int>*)key->ReadObj())[0];
      else if (keyname=="topPtWeightNOM") topPtWeightNOM += (*(vector<double>*)key->ReadObj())[0];
    }
    weight0 *= total / topPtWeightNOM;
  }*/

  //Histograms//
  TString hname;
  hname = Form("costoplep");
  m_Histos1D[hname] = new TH1D(hname,hname,40,-1,1);
  hname = Form("cosantitoplep");
  m_Histos1D[hname] = new TH1D(hname,hname,40,-1,1);
  hname = Form("costopvscosantitop");
  m_Histos2D[hname] = new TH2D(hname,hname,40,-1,1,40,-1,1);
  /*int nDirs = 8;
  for (int i=0; i<nDirs; i++) {
    TString hname = Form("%i_nJet",i);
    m_Histos1D[hname] = new TH1D(hname,hname,MAXJET,0,MAXJET);
    hname = Form("%i_nGoodJet",i);
    m_Histos1D[hname] = new TH1D(hname,hname,MAXJET,0,MAXJET);
    hname = Form("%i_nJetDiff",i);
    m_Histos1D[hname] = new TH1D(hname,hname,MAXJET,0,MAXJET);
    hname = Form("%i_jet0pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
    hname = Form("%i_jet1pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
    hname = Form("%i_jet0eta",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
    hname = Form("%i_jet1eta",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
    hname = Form("%i_jet0phi",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-4,4);
    hname = Form("%i_jet1phi",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-4,4);
    hname = Form("%i_jet0btag",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1);
    hname = Form("%i_jet1btag",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1);
    hname = Form("%i_nbtag",i);
    m_Histos1D[hname] = new TH1D(hname,hname,5,0,5);
    hname = Form("%i_jethT",i);
    m_Histos1D[hname] = new TH1D(hname,hname,300,0,3000);
    hname = Form("%i_jetflavor",i);
    m_Histos1D[hname] = new TH1D(hname,hname,50,-25,25);

    hname = Form("%i_minjet0pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
    hname = Form("%i_minjet1pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
    hname = Form("%i_cleanjet0pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
    hname = Form("%i_cleanjet1pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
    hname = Form("%i_masslmin0",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1000);
    hname = Form("%i_masslmin1",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1000);
    hname = Form("%i_masslljjm",i);
    m_Histos1D[hname] = new TH1D(hname,hname,500,0,5000);

    hname = Form("%i_nEle",i);
    m_Histos1D[hname] = new TH1D(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nEleDiff",i);
    m_Histos1D[hname] = new TH1D(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nMuon",i);
    m_Histos1D[hname] = new TH1D(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nMuonDiff",i);
    m_Histos1D[hname] = new TH1D(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nGoodEle",i);
    m_Histos1D[hname] = new TH1D(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_nGoodMuon",i);
    m_Histos1D[hname] = new TH1D(hname,hname,MAXLEP,0,MAXLEP);
    hname = Form("%i_lep0pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,500,0,5000);
    hname = Form("%i_lep0eta",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
    hname = Form("%i_lep1pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,500,0,5000);
    hname = Form("%i_lep1eta",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
    hname = Form("%i_lep0phi",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-4,4);
    hname = Form("%i_lep1phi",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-4,4);
    hname = Form("%i_dilepmass",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
    hname = Form("%i_lepept",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1000);
    hname = Form("%i_lepmpt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1000);
    hname = Form("%i_deta_lep",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
    hname = Form("%i_deta_lepJet",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
    hname = Form("%i_dphi_jet0met",i);
    m_Histos1D[hname] = new TH1D(hname,hname,70,0,3.5);
    hname = Form("%i_dphi_jet1met",i);
    m_Histos1D[hname] = new TH1D(hname,hname,70,0,3.5);

    hname = Form("%i_muonD0",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-1,1);
    hname = Form("%i_muonDz",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,-1,1);
    hname = Form("%i_rmin0",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,5);
    hname = Form("%i_rmin1",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,5);
    hname = Form("%i_sumrmin",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,10);
    hname = Form("%i_rbl",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,5);
    hname = Form("%i_rl0l1",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,5);
    hname = Form("%i_rl0cleanj",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,5);
    hname = Form("%i_rl1cleanj",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,5);
    hname = Form("%i_lep0perp",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,200);
    hname = Form("%i_lep1perp",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,200);
    hname = Form("%i_lep0perp_in",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,200);
    hname = Form("%i_lep1perp_in",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,200);

    hname = Form("%i_metpt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
    hname = Form("%i_metcorrpt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
    hname = Form("%i_sT",i);
    m_Histos1D[hname] = new TH1D(hname,hname,500,0,5000);
    hname = Form("%i_sT_met",i);
    m_Histos1D[hname] = new TH1D(hname,hname,500,0,5000);
    hname = Form("%i_nPV",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,100);
  }

  TString hname = "muTrigSf";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0.5,1.5);
  hname = "muIdSf";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,2);
  hname = "muTrackSf";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,2);
  hname = "eTrigSf";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,2);
  hname = "eRecoSf";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,2);
  hname = "eIdSf";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,2);

  hname = "jetPt_b";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  hname = "jetPt_c";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  hname = "jetPt_udsg";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  hname = "jetPt_bTagL_b";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  hname = "jetPt_bTagL_c";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  hname = "jetPt_bTagL_udsg";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  hname = "jetPt_bTagM_b";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  hname = "jetPt_bTagM_c";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  hname = "jetPt_bTagM_udsg";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);

  hname = "jetEta_b";
  m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
  hname = "jetEta_c";
  m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
  hname = "jetEta_udsg";
  m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
  hname = "jetEta_bTagL_b";
  m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
  hname = "jetEta_bTagL_c";
  m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
  hname = "jetEta_bTagL_udsg";
  m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
  hname = "jetEta_bTagM_b";
  m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
  hname = "jetEta_bTagM_c";
  m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);
  hname = "jetEta_bTagM_udsg";
  m_Histos1D[hname] = new TH1D(hname,hname,100,-5,5);

  hname = "jetPtDR_b";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,5,200,0,2000);
  hname = "jetPtDR_c";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,5,200,0,2000);
  hname = "jetPtDR_udsg";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,5,200,0,2000);
  hname = "jetPtDR_bTagL_b";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,5,200,0,2000);
  hname = "jetPtDR_bTagL_c";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,5,200,0,2000);
  hname = "jetPtDR_bTagL_udsg";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,5,200,0,2000);
  hname = "jetPtDR_bTagM_b";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,5,200,0,2000);
  hname = "jetPtDR_bTagM_c";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,5,200,0,2000);
  hname = "jetPtDR_bTagM_udsg";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,5,200,0,2000);

  hname = "t_pt";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  hname = "tbar_pt";
  m_Histos1D[hname] = new TH1D(hname,hname,200,0,2000);
  */
  //Set Branches//

  ULong64_t event;
  int run, nPV;
  float rho, mu;

  T->SetBranchAddress("event", &event);
  T->SetBranchAddress("run", &run);
  T->SetBranchAddress("nPV", &nPV);
  T->SetBranchAddress("rho", &rho);
  if (isMC) T->SetBranchAddress("mu", &mu);

 /* string triggers[nTriggers] = {
    "HLT_Mu45_eta2p1_v",
    "HLT_Mu50_v",
    "HLT_TkMu50_v",
    "HLT_Mu30_TkMu11_v",
    "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",
    "HLT_Ele105_CaloIdVT_GsfTrkIdT_v",
    "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_"
  };*/

  vector<bool> *trig_passed = 0;
  //vector<int> *trig_prescale = 0;
  //vector<string> *trig_name = 0;

  T->SetBranchAddress("trig_passed", &trig_passed);
  //T->SetBranchAddress("trig_prescale", &trig_prescale);
  //T->SetBranchAddress("trig_name", &trig_name);

  vector<float> *wgt_env = 0, *wgt_rep = 0;
  int nGen=MAXGEN, gen_status[nGen], gen_PID[nGen], gen_index[nGen], gen_mother0[nGen], gen_mother1[nGen];
  float gen_pt[nGen], gen_mass[nGen], gen_eta[nGen], gen_phi[nGen];

  if (inName.Contains("ttbar", TString::kIgnoreCase)) { //|| inName.Contains("zprime", TString::kIgnoreCase) || inName.Contains("gluon", TString::kIgnoreCase)) {
    T->SetBranchAddress("wgt_env", &wgt_env);
    T->SetBranchAddress("wgt_rep", &wgt_rep);
  }
 if (name.Contains("tt", TString::kIgnoreCase) || name.Contains("zprime", TString::kIgnoreCase) || name.Contains("gkk", TString::kIgnoreCase)) {
    T->SetBranchAddress("nGen", &nGen);
    T->SetBranchAddress("gen_status", gen_status);
    T->SetBranchAddress("gen_PID", gen_PID);
    T->SetBranchAddress("gen_pt", gen_pt);
    T->SetBranchAddress("gen_mass", gen_mass);
    T->SetBranchAddress("gen_eta", gen_eta);
    T->SetBranchAddress("gen_phi", gen_phi);
    T->SetBranchAddress("gen_index",  gen_index);
    T->SetBranchAddress("gen_mother0",  gen_mother0);
    T->SetBranchAddress("gen_mother1",  gen_mother1);
  }

  char lep0flavor, lep1flavor;
  T->SetBranchAddress("lep0flavor", &lep0flavor);
  T->SetBranchAddress("lep1flavor", &lep1flavor);

  int nMuon=MAXLEP, muon_charge[nMuon];
  float muon_eta[nMuon], muon_pt[nMuon], muon_phi[nMuon];
  bool muon_IsTightID[nMuon];

  T->SetBranchAddress("nMuon", &nMuon);
  T->SetBranchAddress("muon_charge", muon_charge);
  T->SetBranchAddress("muon_eta", muon_eta);
  T->SetBranchAddress("muon_pt", muon_pt);
  T->SetBranchAddress("muon_phi", muon_phi);
  T->SetBranchAddress("muon_IsTightID", muon_IsTightID);

  int nEle=MAXLEP, ele_charge[nEle];
  float ele_eta[nEle], ele_pt[nEle], ele_phi[nEle], ele_etaSupClust[nEle];
  bool ele_TightID[nEle];

  T->SetBranchAddress("nEle", &nEle);
  T->SetBranchAddress("ele_charge", ele_charge);
  T->SetBranchAddress("ele_eta", ele_eta);
  T->SetBranchAddress("ele_pt", ele_pt);
  T->SetBranchAddress("ele_phi", ele_phi);
  T->SetBranchAddress("ele_etaSupClust", ele_etaSupClust);
  T->SetBranchAddress("ele_TightID", ele_TightID);

  int nJet=MAXJET;
  int jet_hadflavor[nJet];
  float jet_eta[nJet], jet_phi[nJet], jet_pt[nJet], jet_mass[nJet], jet_area[nJet];
  float jet_btag[nJet], jet_nhf[nJet], jet_nef[nJet], jet_chf[nJet], jet_muf[nJet], jet_elef[nJet], jet_numneutral[nJet], jet_chmult[nJet];
  char jet_clean[nJet];

  T->SetBranchAddress("nJet", &nJet);
  T->SetBranchAddress("jet_eta", jet_eta);
  T->SetBranchAddress("jet_phi", jet_phi);
  T->SetBranchAddress("jet_pt", jet_pt);
  T->SetBranchAddress("jet_mass", jet_mass);
  T->SetBranchAddress("jet_area", jet_area);
  T->SetBranchAddress("jet_clean", jet_clean);

  T->SetBranchAddress("jet_btag", jet_btag);
  T->SetBranchAddress("jet_nhf", jet_nhf);
  T->SetBranchAddress("jet_nef", jet_nef);
  T->SetBranchAddress("jet_chf", jet_chf);
  T->SetBranchAddress("jet_muf", jet_muf);
  T->SetBranchAddress("jet_elef", jet_elef);
  T->SetBranchAddress("jet_numneutral", jet_numneutral);
  T->SetBranchAddress("jet_chmult", jet_chmult);

  int nGenJet=MAXJET;
  float genJet_pt[nGenJet], genJet_eta[nGenJet], genJet_phi[nGenJet], genJet_mass[nGenJet];

  if (isMC) {
    T->SetBranchAddress("jet_hadflavor", jet_hadflavor);

    T->SetBranchAddress("nGenJet", &nGenJet);
    T->SetBranchAddress("genJet_pt", genJet_pt);
    T->SetBranchAddress("genJet_eta", genJet_eta);
    T->SetBranchAddress("genJet_phi", genJet_phi);
    T->SetBranchAddress("genJet_mass", genJet_mass);
  }

  float met_pt, met_px, met_py, met_phi;
  T->SetBranchAddress("met_pt", &met_pt);
  T->SetBranchAddress("met_px", &met_px);
  T->SetBranchAddress("met_py", &met_py);
  T->SetBranchAddress("met_phi", &met_phi);

  //Loop Over Entries//
  int sameRlepjet=0;
  time_t start = time(NULL);

  for (Long64_t n=0; n<nEntries; n++) {
    T->GetEntry(n);

    TLorentzVector lep0, lep1;
    weight = weight0;

    if (isMC) {

      TLorentzVector top, antitop;
      TLorentzVector LEPfromTOP, LEPfromANTITOP;
      TVector3 beta_top , beta_antitop;
      if (name.Contains("tt", TString::kIgnoreCase) || name.Contains("zprime", TString::kIgnoreCase) || name.Contains("gkk", TString::kIgnoreCase)) {
        for (int i=0; i<nGen; i++) {
          if (gen_PID[i]==6 && gen_status[i]>30) top.SetPtEtaPhiM( gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i] );
          else if (gen_PID[i]==-11 || gen_PID[i]==-13) LEPfromTOP.SetPtEtaPhiM( gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i] );
          else if (gen_PID[i]==-6 && gen_status[i]>30) antitop.SetPtEtaPhiM( gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i] );
          else if (gen_PID[i]==11 || gen_PID[i]==13)   LEPfromANTITOP.SetPtEtaPhiM( gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i] );
          beta_top = top.BoostVector() ;
          beta_antitop = antitop.BoostVector() ;
        }
      } 
      double top_lep_cos, antitop_lep_cos;
      if (top.Mag()!=0 && LEPfromTOP.Mag()!=0){
        LEPfromTOP.Boost(-beta_top);
        top_lep_cos =  LEPfromTOP.Vect().Unit().Dot(-beta_top.Unit());
      }
      if (antitop.Mag()!=0 && LEPfromANTITOP.Mag()!=0){
        LEPfromANTITOP.Boost(-beta_antitop);
        antitop_lep_cos =  LEPfromANTITOP.Vect().Unit().Dot(-beta_antitop.Unit());
      }
      FillHist1D("costoplep", top_lep_cos, weight);
      FillHist1D("cosantitoplep", antitop_lep_cos, weight);
      FillHist2D("costopvscosantitop", top_lep_cos, antitop_lep_cos, 1.);
      /*weight *= pileup_weights->GetBinContent( pileup_weights->FindBin(mu) );

      //topPt, pdf, and q2 reweighting - only affect ttbar
      if ( inName.Contains("ttbar", TString::kIgnoreCase) ) {

        //topPt reweighting - use last t's
        //could store this value in branch using ttbar.c to save time in analyze
        if (topPtWeight == "NOMINAL") {
          double t_pt2=0, tbar_pt2=0;
          for (int i=0; i<nGen; i++) {
            if (gen_PID[i]==6 && gen_status[i]>30) t_pt2 = gen_pt[i];
            else if (gen_PID[i]==-6 && gen_status[i]>30) tbar_pt2 = gen_pt[i];
          }
          weight *= sqrt( exp(0.0615-0.0005*t_pt2) * exp(0.0615-0.0005*tbar_pt2) );

          FillHist1D("t_pt", t_pt2, weight);
          FillHist1D("tbar_pt", tbar_pt2, weight);
        }
        //pdf reweighting
        if (pdf != "NOMINAL") {
          vector<float> pdf_plus, pdf_minus;

          for (unsigned int i=0, n=wgt_rep->size(); i<n; i++) {
            if (wgt_rep->at(i) >= 1.) pdf_plus.push_back(wgt_rep->at(i));
            else                      pdf_minus.push_back(wgt_rep->at(i));
          }
          if      (pdf == "UP")   weight *= (1. + rms_pm(pdf_plus));
          else if (pdf == "DOWN") weight *= (1. - rms_pm(pdf_minus));
        }
        //q2 scale reweighting
        if      (q2ttbar == "UP")   weight *= TMath::MaxElement(wgt_env->size(), &wgt_env->at(0));
        else if (q2ttbar == "DOWN") weight *= TMath::MinElement(wgt_env->size(), &wgt_env->at(0));
      }
/*
      //leptonic, semi-leptonic, or hadronic channel from gen information
      if (inName.Contains("ttbar", TString::kIgnoreCase) || inName.Contains("gluon", TString::kIgnoreCase) || inName.Contains("zprime", TString::kIgnoreCase)) {
        TString w_indices = "";
        for (int i=0; i<nGen; i++) { if (abs(gen_PID[i]) == 24) w_indices += to_string(gen_index[i]) + ", "; }
        int num_wl=0;
        for (int i=0; i<nGen; i++) {
          if ( w_indices.Contains( to_string(gen_mother0[i]) + "," ) &&
             ( (abs(gen_PID[i]) == 11 || abs(gen_PID[i]) == 13 || abs(gen_PID[i]) == 15) ) ) num_wl++;
        }
        if      (num_wl == 2) v_cuts_ptr = &v_llcuts;
        else if (num_wl == 1) v_cuts_ptr = &v_ljcuts;
        else                  v_cuts_ptr = &v_jjcuts;
      }
*/
    }
    /*if (channel == "mm") {
      if (lep0flavor == 'm' && lep1flavor == 'm') {
        v_cuts[channelCut].second += weight;  v_cuts_ptr->at(channelCut).second += weight;

        //HLT_Mu50 or HLT_TkMu50 triggers
        if ( !(*trig_passed)[1] && !(*trig_passed)[2] ) continue;

        if (isMC) {
          double muTrackSf = muTrackSfGraph->Eval(muon_eta[0]) * muTrackSfGraph->Eval(muon_eta[1]);
          weight *= muTrackSf;

          int bin0 = muTrigSfHist->FindBin( fabs(muon_eta[0]), muon_pt[0]>muTrig_pT?muTrig_pT:muon_pt[0] );
          //int bin1 = muTrigSfHist->FindBin( fabs(muon_eta[1]), muon_pt[1]>muTrig_pT?muTrig_pT:muon_pt[1] );

          double muTrigSf0 = muTrigSfHist->GetBinContent( bin0 );
          //double muTrigSf1 = muTrigSfHist->GetBinContent( bin1 );
          if      (muTrigSys == "UP")   {
            muTrigSf0 += sqrt( muTrigSfHist->GetBinError( bin0 )*muTrigSfHist->GetBinError( bin0 ) + 0.01*muTrigSf0*0.01*muTrigSf0 );
            //muTrigSf1 += sqrt( muTrigSfHist->GetBinError( bin1 )*muTrigSfHist->GetBinError( bin1 ) + 0.01*muTrigSf1*0.01*muTrigSf1 );
          }
          else if (muTrigSys == "DOWN") {
            muTrigSf0 -= sqrt( muTrigSfHist->GetBinError( bin0 )*muTrigSfHist->GetBinError( bin0 ) + 0.01*muTrigSf0*0.01*muTrigSf0 );
            //muTrigSf1 -= sqrt( muTrigSfHist->GetBinError( bin1 )*muTrigSfHist->GetBinError( bin1 ) + 0.01*muTrigSf1*0.01*muTrigSf1 );
          }
          if (muon_pt[0] < 53) muTrigSf0=1;
          //if (muon_pt[1] < 53) muTrigSf1=1;
          //weight *= muTrigSf0*muTrigSf1;
          weight *= muTrigSf0;

          FillHist1D("muTrackSf", muTrackSf, 1.);
          FillHist1D("muTrigSf", muTrigSf0, 1.);
        }
        v_cuts[trigCut].second += weight;  v_cuts_ptr->at(trigCut).second += weight;

        if ( !muon_IsTightID[0] || !muon_IsTightID[1] ) continue;
        if ( muon_pt[0] < 53 || muon_pt[1] < 25 ) continue;
        if (fabs(muon_eta[0]) > 2.4 || fabs(muon_eta[1]) > 2.4) continue;

        if (isMC) {
          int bin0 = muIdSfHist->FindBin( fabs(muon_eta[0]), muon_pt[0]>muId_pT?muId_pT:muon_pt[0] );
          int bin1 = muIdSfHist->FindBin( fabs(muon_eta[1]), muon_pt[1]>muId_pT?muId_pT:muon_pt[1] );

          double muIdSf0 = muIdSfHist->GetBinContent( bin0 );
          double muIdSf1 = muIdSfHist->GetBinContent( bin1 );
          if      (muIdSys == "UP")   {
            muIdSf0 += sqrt( muIdSfHist->GetBinError( bin0 )*muIdSfHist->GetBinError( bin0 ) + 0.02*muIdSf0*0.02*muIdSf0 );
            muIdSf1 += sqrt( muIdSfHist->GetBinError( bin1 )*muIdSfHist->GetBinError( bin1 ) + 0.02*muIdSf1*0.02*muIdSf1 );
          }
          else if (muIdSys == "DOWN") {
            muIdSf0 -= sqrt( muIdSfHist->GetBinError( bin0 )*muIdSfHist->GetBinError( bin0 ) + 0.02*muIdSf0*0.02*muIdSf0 );
            muIdSf1 -= sqrt( muIdSfHist->GetBinError( bin1 )*muIdSfHist->GetBinError( bin1 ) + 0.02*muIdSf1*0.02*muIdSf1 );
          }

          weight *= muIdSf0*muIdSf1;
          FillHist1D("muIdSf", muIdSf0*muIdSf1, 1.);
        }
        v_cuts[lepkinCut].second += weight;  v_cuts_ptr->at(lepkinCut).second += weight;

        if (muon_charge[0]*muon_charge[1] > 0) continue;
        v_cuts[signCut].second += weight;  v_cuts_ptr->at(signCut).second += weight;

        //use these events for em channel
        if ( nEle>0 && ele_TightID[0] && ele_pt[0]>25 && fabs(ele_eta[0])<2.5 ) continue;
        v_cuts[thirdLepCut].second += weight;  v_cuts_ptr->at(thirdLepCut).second += weight;

        lep0.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
        lep1.SetPtEtaPhiM(muon_pt[1], muon_eta[1], muon_phi[1], MUONMASS);

        if ((lep0+lep1).M() < 20) continue;
      }
      else continue;
    }
    else if (channel == "ee") {
      if (lep0flavor == 'e' && lep1flavor == 'e') {
        v_cuts[channelCut].second += weight;  v_cuts_ptr->at(channelCut).second += weight;

        //HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_ trigger
        if ( !(*trig_passed)[7] ) continue;

        if (isMC) {
          int bin0 = eTrigSfHist->FindBin( fabs(ele_etaSupClust[0]), ele_pt[0]>eTrig_pT?eTrig_pT-1:ele_pt[0] );
          int bin1 = eTrigSfHist->FindBin( fabs(ele_etaSupClust[1]), ele_pt[1]>eTrig_pT?eTrig_pT-1:ele_pt[1] );

          double eTrigSf0 = eTrigSfHist->GetBinContent( bin0 );
          double eTrigSf1 = eTrigSfHist->GetBinContent( bin1 );
          if      (eleTrigSys == "UP")   {
            double af = (ele_pt[0]>eTrig_pT) ? 2 : 1;
            eTrigSf0 += sqrt( af*eTrigSfHist->GetBinError( bin0 )*af*eTrigSfHist->GetBinError( bin0 ) + 0.02*eTrigSf0*0.02*eTrigSf0 );

            af = (ele_pt[1]>eTrig_pT) ? 2 : 1;
            eTrigSf1 += af*eTrigSfHist->GetBinError( bin1 );
          }
          else if (eleTrigSys == "DOWN") {
            double af = (ele_pt[0]>eTrig_pT) ? 2 : 1;
            eTrigSf0 -= sqrt( af*eTrigSfHist->GetBinError( bin0 )*af*eTrigSfHist->GetBinError( bin0 ) + 0.02*eTrigSf0*0.02*eTrigSf0 );

            af = (ele_pt[1]>eTrig_pT) ? 2 : 1;
            eTrigSf1 -= af*eTrigSfHist->GetBinError( bin1 );
          }
          weight *= eTrigSf0*eTrigSf1;

          double eRecoSf = eRecoSfHist->GetBinContent( eRecoSfHist->FindBin( ele_etaSupClust[0], ele_pt[0]>eReco_pT?eReco_pT:ele_pt[0] ) )
                         * eRecoSfHist->GetBinContent( eRecoSfHist->FindBin( ele_etaSupClust[1], ele_pt[1]>eReco_pT?eReco_pT:ele_pt[1] ) );
          weight *= eRecoSf;

          FillHist1D("eTrigSf", eTrigSf0*eTrigSf1, 1.);
          FillHist1D("eRecoSf", eRecoSf, 1.);
        }
        v_cuts[trigCut].second += weight;  v_cuts_ptr->at(trigCut).second += weight;

        if ( !ele_TightID[0] || !ele_TightID[1] ) continue;
        if ( ele_pt[0] < 45 || ele_pt[1] < 36 ) continue;
        if (fabs(ele_eta[0]) > 2.5 || fabs(ele_eta[1]) > 2.5) continue;

        if (isMC) {
          int bin0 = eIdSfHist->FindBin( ele_etaSupClust[0], ele_pt[0]>eId_pT?eId_pT:ele_pt[0] );
          int bin1 = eIdSfHist->FindBin( ele_etaSupClust[1], ele_pt[1]>eId_pT?eId_pT:ele_pt[1] );

          double eIdSf0 = eIdSfHist->GetBinContent( bin0 );
          double eIdSf1 = eIdSfHist->GetBinContent( bin1 );
          if      (eleIdSys == "UP")   {
            eIdSf0 += eIdSfHist->GetBinError( bin0 );
            eIdSf1 += eIdSfHist->GetBinError( bin1 );
          }
          else if (eleIdSys == "DOWN") {
            eIdSf0 -= eIdSfHist->GetBinError( bin0 );
            eIdSf1 -= eIdSfHist->GetBinError( bin1 );
          }

          weight *= eIdSf0*eIdSf1;
          FillHist1D("eIdSf",  eIdSf0*eIdSf1, 1.);
        }
        v_cuts[lepkinCut].second += weight;  v_cuts_ptr->at(lepkinCut).second += weight;

        if (ele_charge[0]*ele_charge[1] > 0) continue;
        v_cuts[signCut].second += weight;  v_cuts_ptr->at(signCut).second += weight;

        //use these events for em channel
        if ( nMuon>0 && muon_IsTightID[0] && muon_pt[0]>53 && fabs(muon_eta[0])<2.4 ) continue;
        v_cuts[thirdLepCut].second += weight;  v_cuts_ptr->at(thirdLepCut).second += weight;

        lep0.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
        lep1.SetPtEtaPhiM(ele_pt[1], ele_eta[1], ele_phi[1], ELEMASS);

        if ((lep0+lep1).M() < 20) continue;
      }
      else continue;
    }
    else {
      if (lep0flavor != lep1flavor) {
        v_cuts[channelCut].second += weight;  v_cuts_ptr->at(channelCut).second += weight;

        //HLT_Mu50 or HLT_TkMu50 triggers
        if ( !(*trig_passed)[1] && !(*trig_passed)[2] ) continue;

        if (isMC) {
          double muTrackSf = muTrackSfGraph->Eval(muon_eta[0]);

          int bin = muTrigSfHist->FindBin( fabs(muon_eta[0]), muon_pt[0]>muTrig_pT?muTrig_pT:muon_pt[0] );
          double muTrigSf = muTrigSfHist->GetBinContent( bin );
          if      (muTrigSys == "UP")   muTrigSf += sqrt( muTrigSfHist->GetBinError( bin )*muTrigSfHist->GetBinError( bin ) + 0.01*muTrigSf*0.01*muTrigSf );
          else if (muTrigSys == "DOWN") muTrigSf -= sqrt( muTrigSfHist->GetBinError( bin )*muTrigSfHist->GetBinError( bin ) + 0.01*muTrigSf*0.01*muTrigSf );
          if (muon_pt[0] < 53) muTrigSf=1;

          weight *= muTrackSf * muTrigSf;

          FillHist1D("muTrigSf", muTrigSf, 1.);
          FillHist1D("muTrackSf", muTrackSf, 1.);
        }
        v_cuts[trigCut].second += weight;  v_cuts_ptr->at(trigCut).second += weight;

        if ( !muon_IsTightID[0] || !ele_TightID[0] ) continue;
        if ( muon_pt[0] < 53 || ele_pt[0] < 25 ) continue;
        if (fabs(muon_eta[0]) > 2.4 || fabs(ele_eta[0]) > 2.5) continue;

        if (isMC) {
          int bin = muIdSfHist->FindBin( fabs(muon_eta[0]), muon_pt[0]>muId_pT?muId_pT:muon_pt[0] );
          double muIdSf = muIdSfHist->GetBinContent( bin );
          if      (muIdSys == "UP")   muIdSf += sqrt( muIdSfHist->GetBinError( bin )*muIdSfHist->GetBinError( bin ) + 0.02*muIdSf*0.02*muIdSf );
          else if (muIdSys == "DOWN") muIdSf -= sqrt( muIdSfHist->GetBinError( bin )*muIdSfHist->GetBinError( bin ) + 0.02*muIdSf*0.02*muIdSf );

          double eRecoSf = eRecoSfHist->GetBinContent( eRecoSfHist->FindBin( ele_etaSupClust[0], ele_pt[0]>eReco_pT?eReco_pT:ele_pt[0] ) );

          bin = eIdSfHist->FindBin( ele_etaSupClust[0], ele_pt[0]>eId_pT?eId_pT:ele_pt[0] );
          double eIdSf = eIdSfHist->GetBinContent( bin );
          if      (eleIdSys == "UP")   eIdSf += eIdSfHist->GetBinError( bin );
          else if (eleIdSys == "DOWN") eIdSf -= eIdSfHist->GetBinError( bin );

          weight *= muIdSf * eRecoSf * eIdSf;

          FillHist1D("muIdSf", muIdSf, 1.);
          FillHist1D("eRecoSf", eRecoSf, 1.);
          FillHist1D("eIdSf", eIdSf, 1.);
        }

        v_cuts[lepkinCut].second += weight;  v_cuts_ptr->at(lepkinCut).second += weight;

        if (ele_charge[0]*muon_charge[0] > 0) continue;
        v_cuts[signCut].second += weight;  v_cuts_ptr->at(signCut).second += weight;

        v_cuts[thirdLepCut].second += weight;  v_cuts_ptr->at(thirdLepCut).second += weight;

        if (lep0flavor=='e') {
          lep0.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
          lep1.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
        }
        else {
          lep0.SetPtEtaPhiM(muon_pt[0], muon_eta[0], muon_phi[0], MUONMASS);
          lep1.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ELEMASS);
        }
      }
      else continue;
    }
    v_cuts[dilepmassCut].second += weight;  v_cuts_ptr->at(dilepmassCut).second += weight;
    double dilepmass = (lep0+lep1).M();

    if ( lep0flavor==lep1flavor && (76<dilepmass && dilepmass<106) ) continue;
    v_cuts[dilepVetoCut].second += weight;  v_cuts_ptr->at(dilepVetoCut).second += weight;

    string era = eras[0];
    for( map<string, pair<int, int> >::const_iterator it = m_IOV.begin(); it != m_IOV.end(); ++it ) {
      const pair<int, int>& interval = it->second;
      if (interval.first <= run && run <= interval.second) { era = it->first; break; }
    }

    vector<pair<int, double> > jet_index_corrpt;
    TLorentzVector minjet0, minjet1;
    double rmin0=99, rmin1=99;
    int minjet0idx=-1, minjet1idx=-1;
    double ctype1_x=0, ctype1_y=0;
    double rl0cleanj=-1, rl1cleanj=-1, cleanjet0pt=-1, cleanjet1pt=-1;

    TRandom3* rand = new TRandom3(event);
    vector<float> rands (nJet, 0.);
    double hT=0;
    for (int i=0; i<nJet; i++) {

      //loose jet cut
      if (fabs(jet_eta[i]) <= 2.7) {
        if (jet_nhf[i]>=0.99 || jet_nef[i]>=0.99 || (jet_numneutral[i]+jet_chmult[i])<=1) continue;
        if (fabs(jet_eta[i]) <= 2.4 && ( jet_chf[i]<=0 || jet_chmult[i]<=0 || jet_elef[i]>=0.99 )) continue;
      }
      else if (2.7 < fabs(jet_eta[i]) && fabs(jet_eta[i]) <= 3.0) {
        if (jet_nhf[i]>=0.98 || jet_nef[i]<=0.01 || jet_numneutral[i]<=2) continue;
      }
      else {
        if (jet_nef[i]>=0.9 || jet_numneutral[i]<=10) continue;
      }

      jetCorrectors[era]->setJetEta( jet_eta[i] );
      jetCorrectors[era]->setJetPt( jet_pt[i] );
      jetCorrectors[era]->setJetA( jet_area[i] );
      jetCorrectors[era]->setRho(rho);

      TLorentzVector jet;
      jet.SetPtEtaPhiM(jet_pt[i], jet_eta[i], jet_phi[i], jet_mass[i]);
      jet *= jetCorrectors[era]->getCorrection();

      if (isMC) {
        rands[i] = rand->Uniform(1.);

        if (jec=="UP" || jec=="DOWN") {
          double sign = jec=="UP" ? +1. : -1 ;
          jecUncert[era]->setJetEta( jet.Eta() );
          jecUncert[era]->setJetPt( jet.Pt() );
          double unc = jecUncert[era]->getUncertainty(true);
          jet *= 1. + sign*unc ;
        }

        JME::JetParameters res_pars;
        res_pars.setJetEta( jet.Eta() );
        res_pars.setJetPt( jet.Pt() ); //corrected pt
        res_pars.setRho(rho);

        double jet_res = res_obj.getResolution( res_pars );
        Variation jer_sys = Variation::NOMINAL;
        if (jer=="UP") jer_sys = Variation::UP;
        else if (jer=="DOWN") jer_sys = Variation::DOWN;

        double jet_ressf = ressf_obj.getScaleFactor( res_pars , jer_sys );

        TLorentzVector matched_genJet;
        double rmin_genJet = 99.;
        for (int i_gen=0; i_gen<nGenJet; i_gen++) {

          TLorentzVector genJet;
          genJet.SetPtEtaPhiM(genJet_pt[i_gen], genJet_eta[i_gen], genJet_phi[i_gen], genJet_mass[i_gen]);

          double dR = jet.DeltaR(genJet);
          if ( dR<rmin_genJet && dR<0.2 && fabs(jet.Pt()-genJet.Pt()) <= 3*jet_res*jet.Pt() ) {

            rmin_genJet = dR;
            matched_genJet = genJet;
          }
        }

        double smearFactor = 1.;
        if (matched_genJet.E() != 0) smearFactor = 1. + (jet_ressf - 1.) * (jet.Pt() - matched_genJet.Pt()) / jet.Pt();
        else if (jet_ressf > 1) {
          double sigma = jet_res * sqrt(jet_ressf * jet_ressf - 1.);
          smearFactor = 1. + rand->Gaus(0, sigma);
        }
        if (smearFactor > 0) jet *= smearFactor;
      }

      //corrected MET
      if ( jet.Pt()>15 && (jet_elef[i]+jet_nef[i])<0.9 ) {
        jetL1Correctors[era]->setJetEta( jet_eta[i] );
        jetL1Correctors[era]->setJetPt( jet_pt[i] ); //uncorrected pt
        jetL1Correctors[era]->setJetA( jet_area[i] );
        jetL1Correctors[era]->setRho(rho);

        TLorentzVector jetL1;
        jetL1.SetPtEtaPhiM( jet_pt[i], jet_eta[i], jet_phi[i], jet_mass[i] );
        jetL1 *= jetL1Correctors[era]->getCorrection();

        ctype1_x += (jet.Px()-jetL1.Px());
        ctype1_y += (jet.Py()-jetL1.Py());
      }

      if (jet.Pt()>15 && fabs(jet_eta[i])<3.) {

        if (lep0.DeltaR(jet) < rmin0) {
          rmin0 = lep0.DeltaR(jet);
          minjet0 = jet;
          minjet0idx = i;
        }
        if (lep1.DeltaR(jet) < rmin1) {
          rmin1 = lep1.DeltaR(jet);
          minjet1 = jet;
          minjet1idx = i;
        }
        if (jet_clean[i] == 'l' || jet_clean[i] == 'b') { rl0cleanj = lep0.DeltaR(jet); cleanjet0pt = jet.Pt(); }
        if (jet_clean[i] == 's' || jet_clean[i] == 'b') { rl1cleanj = lep1.DeltaR(jet); cleanjet1pt = jet.Pt(); }

        if (jet.Pt()>30 && fabs(jet_eta[i])<2.4) {
          jet_index_corrpt.push_back( make_pair(i, jet.Pt()) );
          hT+=jet.Pt();
        }
      }
    }
    int nGoodJet = jet_index_corrpt.size();
    if (nGoodJet < 2) continue;

    if (minjet0 == minjet1) sameRlepjet++;
    double lep0perp = lep0.Perp( minjet0.Vect() );
    double lep1perp = lep1.Perp( minjet1.Vect() );

    if( (lep0perp<15 && rmin0<0.4) || (lep1perp<15 && rmin1<0.4) ) continue;
    v_cuts[ptrelCut].second += weight;  v_cuts_ptr->at(ptrelCut).second += weight;

    if (setSUMDRCut.Contains("ON",TString::kIgnoreCase) ) {
      if (rmin0+rmin1 >= 2.) continue;
      if (setSUMDRCut=="ONbt") { if (rmin0+rmin1 >= 1.) continue; }
      if (setSUMDRCut=="ONnb") { if (rmin0+rmin1 < 1.)  continue; }
    }
    else if (setSUMDRCut=="REVERSE") { if (rmin0+rmin1 < 2.) continue; }
    v_cuts[dRCut].second += weight;  v_cuts_ptr->at(dRCut).second += weight;

    sort(jet_index_corrpt.begin(), jet_index_corrpt.end(), sortJetPt);
    int jet0index = jet_index_corrpt[0].first, jet1index = jet_index_corrpt[1].first;
    double jet0pt = jet_index_corrpt[0].second, jet1pt = jet_index_corrpt[1].second;

    //at least one jet
    if (jet0pt < 100) continue;
    v_cuts[jetCut].second += weight;  v_cuts_ptr->at(jetCut).second += weight;

    double met_corrpx = met_px - ctype1_x;
    double met_corrpy = met_py - ctype1_y;
    double met_corrpt = sqrt(met_corrpx*met_corrpx + met_corrpy*met_corrpy);

    TLorentzVector met;
    met.SetPtEtaPhiE(met_corrpt, 0, met_phi, met_corrpt);
    TLorentzVector jet0, jet1;
    jet0.SetPtEtaPhiM(jet0pt, jet_eta[jet0index], jet_phi[jet0index], jet0pt / jet_pt[jet0index] * jet_mass[jet0index]);
    jet1.SetPtEtaPhiM(jet1pt, jet_eta[jet1index], jet_phi[jet1index], jet1pt / jet_pt[jet1index] * jet_mass[jet1index]);

    int nGoodMuon=0;
    for (int i=0; i<nMuon; i++) { if (muon_IsTightID[i]) nGoodMuon++; }

    int nGoodEle=0;
    for (int i=0; i<nEle; i++) { if (ele_TightID[i]) nGoodEle++; }

    bool jet0btag = jet_btag[jet0index] > btagWP_L;
    bool jet1btag = jet_btag[jet1index] > btagWP_L;

    int jetflavor0=-25, jetflavor1=-25;
    if (isMC) {
      jetflavor0 = jet_hadflavor[jet0index];
      jetflavor1 = jet_hadflavor[jet1index];

      //btag eff for two jets
      if (jet1pt > 50) {

        bool jet0btagM = jet_btag[jet0index] > btagWP_M;
        bool jet1btagM = jet_btag[jet1index] > btagWP_M;

        TString jetflavor0str, jetflavor1str;
        if      ( abs(jetflavor0) == 4 ) jetflavor0str = "c";
        else if ( abs(jetflavor0) == 5 ) jetflavor0str = "b";
        else                             jetflavor0str = "udsg";

        if      ( abs(jetflavor1) == 4 ) jetflavor1str = "c";
        else if ( abs(jetflavor1) == 5 ) jetflavor1str = "b";
        else                             jetflavor1str = "udsg";

        FillHist1D("jetPt_" +jetflavor0str, jet0pt, 1.);
        FillHist1D("jetEta_"+jetflavor0str, jet0.Eta(), 1.);
        FillHist1D("jetPt_" +jetflavor1str, jet1pt, 1.);
        FillHist1D("jetEta_"+jetflavor1str, jet1.Eta(), 1.);

        if (jet0btag) {
          FillHist1D("jetPt_bTagL_" +jetflavor0str, jet0pt, 1.);
          FillHist1D("jetEta_bTagL_"+jetflavor0str, jet0.Eta(), 1.);
        }
        if (jet0btagM) {
          FillHist1D("jetPt_bTagM_" +jetflavor0str, jet0pt, 1.);
          FillHist1D("jetEta_bTagM_"+jetflavor0str, jet0.Eta(), 1.);
        }
        if (jet1btag) {
          FillHist1D("jetPt_bTagL_" +jetflavor1str, jet1pt, 1.);
          FillHist1D("jetEta_bTagL_"+jetflavor1str, jet1.Eta(), 1.);
        }
        if (jet1btagM) {
          FillHist1D("jetPt_bTagM_" +jetflavor1str, jet1pt, 1.);
          FillHist1D("jetEta_bTagM_"+jetflavor1str, jet1.Eta(), 1.);
        }

        double jet0dr=-1, jet1dr=-1;
        if      (jet0index == minjet0idx) jet0dr = rmin0;
        else if (jet0index == minjet1idx) jet0dr = rmin1;
        if      (jet1index == minjet0idx) jet1dr = rmin0;
        else if (jet1index == minjet1idx) jet1dr = rmin1;

        if (jet0dr != -1) {
          FillHist2D("jetPtDR_"+jetflavor0str, jet0dr, jet0pt, 1.);
          if (jet0btag)  FillHist2D("jetPtDR_bTagL_"+jetflavor0str, jet0dr, jet0pt, 1.);
          if (jet0btagM) FillHist2D("jetPtDR_bTagM_"+jetflavor0str, jet0dr, jet0pt, 1.);
        }
        if (jet1dr != -1) {
          FillHist2D("jetPtDR_"+jetflavor1str, jet1dr, jet1pt, 1.);
          if (jet1btag)  FillHist2D("jetPtDR_bTagL_"+jetflavor1str, jet1dr, jet1pt, 1.);
          if (jet1btagM) FillHist2D("jetPtDR_bTagM_"+jetflavor1str, jet1dr, jet1pt, 1.);
        }
      }
      TString variation0 = btagSF, variation1 = btagSF;
      TGraphAsymmErrors* eff0, *eff1;
      if ( abs(jetflavor0) == 4 ) eff0 = btag_eff_c;
      else if ( abs(jetflavor0) == 5 ) eff0 = btag_eff_b;
      else { eff0 = btag_eff_udsg; variation0 = mistagSF; }

      if ( abs(jetflavor1) == 4 ) eff1 = btag_eff_c;
      else if ( abs(jetflavor1) == 5 ) eff1 = btag_eff_b;
      else { eff1 = btag_eff_udsg; variation1 = mistagSF; }

      jet0btag = newBTag( rands[jet0index], jet0pt, jetflavor0, jet0btag, *eff0, variation0 );
      jet1btag = newBTag( rands[jet1index], jet1pt, jetflavor1, jet1btag, *eff1, variation1 );
    }
    delete rand;

    double rl0l1 = lep0.DeltaR(lep1);
    double lepept=0, lepmpt=0;
    if (lep0flavor == 'm') lepmpt += lep0.Pt();
    else lepept += lep0.Pt();
    if (lep1flavor == 'm') lepmpt += lep1.Pt();
    else lepept += lep1.Pt();

    double rbal=-1, rabl=-1;
    if (inName.Contains("ttbar", TString::kIgnoreCase)) {
      //lepton, anti-lepton, b quark, anti-b quark
      TLorentzVector glep, galep, gb, gab;

      for (int i=0; i<nGen; i++) {

        if (gen_PID[i] == 5) gb.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
        else if (gen_PID[i] == -5) gab.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
        else if (gen_PID[i] == 11 || gen_PID[i] == 13) glep.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
        else if (gen_PID[i] == -11 || gen_PID[i] == -13) galep.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
      }

      //b goes with anti-lepton and anti-b goes with lepton
      if (gb.Pt() != 0 && galep.Pt() != 0) rbal = gb.DeltaR(galep);
      if (gab.Pt() != 0 && glep.Pt() != 0) rabl = gab.DeltaR(glep);
    }

    TString prefix;

    //exactly one jet
    if (jet1pt < 50) {
      //zero btags
      if (!jet0btag) {
        if (met_corrpt < 30) { v_cuts[zerobtagCut1jet_metR].second += weight; v_cuts_ptr->at(zerobtagCut1jet_metR).second += weight; prefix="0_"; }
        else                 { v_cuts[zerobtagCut1jet].second += weight;      v_cuts_ptr->at(zerobtagCut1jet).second += weight;      prefix="4_"; }
      }
      //at least one btag
      else {
        if (met_corrpt < 30) { v_cuts[onebtagCut1jet_metR].second += weight; v_cuts_ptr->at(onebtagCut1jet_metR).second += weight; prefix="1_"; }
        else                 { v_cuts[onebtagCut1jet].second += weight;      v_cuts_ptr->at(onebtagCut1jet).second += weight;      prefix="5_"; }
      }
    }
    //at least two jets
    else {
      //exactly zero btags
      if ( !jet0btag && !jet1btag ) {
        if (met_corrpt < 30) { v_cuts[zerobtagCut2jets_metR].second += weight; v_cuts_ptr->at(zerobtagCut2jets_metR).second += weight; prefix="2_"; }
        else                 { v_cuts[zerobtagCut2jets].second += weight;      v_cuts_ptr->at(zerobtagCut2jets).second += weight;      prefix="6_"; }
      }
      //at least one btag
      else {
        if (met_corrpt < 30) { v_cuts[morethan0btagCut2jets_metR].second += weight; v_cuts_ptr->at(morethan0btagCut2jets_metR).second += weight; prefix="3_"; }
        else                 { v_cuts[morethan0btagCut2jets].second += weight;      v_cuts_ptr->at(morethan0btagCut2jets).second += weight;      prefix="7_"; }
      }
    }

    FillHist1D(prefix+"nEleDiff", nEle-nGoodEle, weight);
    FillHist1D(prefix+"nMuonDiff", nMuon-nGoodMuon, weight);
    FillHist1D(prefix+"nJetDiff", nJet-nGoodJet, weight);
    FillHist1D(prefix+"nEle", nEle, weight);
    FillHist1D(prefix+"nMuon", nMuon, weight);
    FillHist1D(prefix+"nJet", nJet, weight);
    FillHist1D(prefix+"nGoodEle", nGoodEle, weight);
    FillHist1D(prefix+"nGoodMuon", nGoodMuon, weight);
    FillHist1D(prefix+"nGoodJet", nGoodJet, weight);

    FillHist1D(prefix+"lep0pt", lep0.Pt(), weight);
    FillHist1D(prefix+"lep0eta", lep0.Eta(), weight);
    FillHist1D(prefix+"lep0phi", lep0.Phi(), weight);
    FillHist1D(prefix+"lep1pt", lep1.Pt(), weight);
    FillHist1D(prefix+"lep1eta", lep1.Eta(), weight);
    FillHist1D(prefix+"lep1phi", lep1.Phi(), weight);
    FillHist1D(prefix+"dilepmass", dilepmass, weight);
    FillHist1D(prefix+"lepept", lepept, weight);
    FillHist1D(prefix+"lepmpt", lepmpt, weight);

    FillHist1D(prefix+"rmin0", rmin0, weight);
    FillHist1D(prefix+"rmin1", rmin1, weight);
    FillHist1D(prefix+"sumrmin", rmin0+rmin1, weight);
    FillHist1D(prefix+"rl0l1", rl0l1, weight);
    FillHist1D(prefix+"rl0cleanj", rl0cleanj, weight);
    FillHist1D(prefix+"rl1cleanj", rl1cleanj, weight);
    FillHist1D(prefix+"lep0perp", lep0perp, weight);
    FillHist1D(prefix+"lep1perp", lep1perp, weight);

    if (rmin0 < 0.4) FillHist1D(prefix+"lep0perp_in", lep0perp, weight);
    if (rmin1 < 0.4) FillHist1D(prefix+"lep1perp_in", lep1perp, weight);

    FillHist1D(prefix+"jet0pt", jet0.Pt(), weight);
    FillHist1D(prefix+"jet0eta", jet0.Eta(), weight);
    FillHist1D(prefix+"jet0phi", jet0.Phi(), weight);
    FillHist1D(prefix+"jet0btag", jet_btag[jet0index], weight);
    FillHist1D(prefix+"jetflavor", jetflavor0, weight);
    FillHist1D(prefix+"jet1pt", jet1.Pt(), weight);
    FillHist1D(prefix+"jet1eta", jet1.Eta(), weight);
    FillHist1D(prefix+"jet1phi", jet1.Phi(), weight);
    FillHist1D(prefix+"jet1btag", jet_btag[jet1index], weight);
    FillHist1D(prefix+"jetflavor", jetflavor1, weight);
    FillHist1D(prefix+"nbtag", int(jet0btag)+int(jet1btag), weight);

    FillHist1D(prefix+"jethT", hT, weight);
    FillHist1D(prefix+"metpt", met_pt, weight);
    FillHist1D(prefix+"metcorrpt", met_corrpt, weight);

    double sT = hT+lep0.Pt()+lep1.Pt();
    double sT_met = sT + met_corrpt;
    double masslljjm = (lep0+lep1+jet0+jet1+met).M();

    FillHist1D(prefix+"sT", sT, weight);
    FillHist1D(prefix+"sT_met", sT_met>=5000 ? 4999.9 : sT_met, weight);
    FillHist1D(prefix+"masslljjm", masslljjm>=5000 ? 4999.9 : masslljjm, weight);

    if (rbal != -1) FillHist1D(prefix+"rbl", rbal, weight);
    if (rabl != -1) FillHist1D(prefix+"rbl", rabl, weight);

    FillHist1D(prefix+"minjet0pt", minjet0.Pt(), weight);
    FillHist1D(prefix+"minjet1pt", minjet1.Pt(), weight);
    FillHist1D(prefix+"cleanjet0pt", cleanjet0pt, weight);
    FillHist1D(prefix+"cleanjet1pt", cleanjet1pt, weight);
    FillHist1D(prefix+"masslmin0", (lep0+minjet0).M(), weight);
    FillHist1D(prefix+"masslmin1", (lep1+minjet1).M(), weight);

    FillHist1D(prefix+"deta_lep", lep0.Eta() - lep1.Eta(), weight);
    FillHist1D(prefix+"deta_lepJet", (lep0+minjet0).Eta() - (lep1+minjet1).Eta(), weight);
    FillHist1D(prefix+"dphi_jet0met", fabs( deltaPhi( jet0.Phi(), met.Phi() ) ), weight);
    FillHist1D(prefix+"dphi_jet1met", fabs( deltaPhi( jet1.Phi(), met.Phi() ) ), weight);

    FillHist1D(prefix+"nPV", nPV, weight);*/
  }
  cout << difftime(time(NULL), start) << " s" << endl;
  cout << "Min_jet0 = Min_jet1: " << sameRlepjet << endl;

  /*TH1D* cuts = new TH1D("cuts","cuts",numCuts,-0.5,float(numCuts)-0.5);
  TH1D* ll_cuts = new TH1D("ll_cuts","ll_cuts",numCuts,-0.5,float(numCuts)-0.5);
  TH1D* lj_cuts = new TH1D("lj_cuts","lj_cuts",numCuts,-0.5,float(numCuts)-0.5);
  TH1D* jj_cuts = new TH1D("jj_cuts","jj_cuts",numCuts,-0.5,float(numCuts)-0.5);

  //Cutflow Table//
  cout<<"===================================================================================================\n";
  cout<<"                                     Cut Flow Table: " + inName( inName.Last('/')+1, inName.Index('.')-inName.Last('/')-1 ) + "\n";
  cout<<"===================================================================================================\n";

  cout<<      "                               |||          Nevent          |||     Efficiency (Relative Efficiency)";

//  if (inName.Contains("ttbar", TString::kIgnoreCase) || inName.Contains("zprime", TString::kIgnoreCase) || inName.Contains("gluon", TString::kIgnoreCase))
//    cout << "   Nllgen   |||      Nljgen   |||     Njjgen\n";
//  else
    cout << "\n";

  for (int i=0; i<numCuts; i++) {
    if (i == 0)
      cout << Form("%-30s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[0].second);

    else if (i > jetCut)
      cout << Form("%-30s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[jetCut].second);

    else
      cout << Form("%-30s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[i-1].second);

//    if (inName.Contains("ttbar", TString::kIgnoreCase) || inName.Contains("zprime", TString::kIgnoreCase) || inName.Contains("gluon", TString::kIgnoreCase))
//      cout << Form("       |||   %9.1f   |||   %9.1f   |||   %9.1f", v_llcuts[i].second, v_ljcuts[i].second, v_jjcuts[i].second) << endl;
//    else
      cout << endl;

    if (i==countMet || i==jetCut)
      cout << "---------------------------------------------------------------------------------------------------" << endl;

    cuts->SetBinContent(i+1, v_cuts[i].second);        ll_cuts->SetBinContent(i+1, v_llcuts[i].second);      lj_cuts->SetBinContent(i+1, v_ljcuts[i].second);      jj_cuts->SetBinContent(i+1, v_jjcuts[i].second);
    cuts->GetXaxis()->SetBinLabel(i+1,Form("%i",i+1)); ll_cuts->GetXaxis()->SetBinLabel(i+1,Form("%i",i+1)); lj_cuts->GetXaxis()->SetBinLabel(i+1,Form("%i",i+1)); jj_cuts->GetXaxis()->SetBinLabel(i+1,Form("%i",i+1));
  }
  cout << endl;
  */
  //Write Histograms//
  outName = name + "_analyzed.root";
  TFile* outFile = new TFile(outName,"RECREATE");
  outFile->cd();

  //cuts->Write(); ll_cuts->Write(); lj_cuts->Write(); jj_cuts->Write();

  //for (int i=0; i<nDirs; i++) outFile->mkdir( Form("%i/", i) );

  for (map<TString, TH1*>::iterator hid = m_Histos1D.begin(); hid != m_Histos1D.end(); hid++) {
    outFile->cd();
    TString prefix = hid->first(0, 1);

    if ( prefix.IsDigit() ) outFile->cd(outName + ":/" + prefix);

    hid->second->Write();
  }
  for (map<TString, TH2*>::iterator hid = m_Histos2D.begin(); hid != m_Histos2D.end(); hid++) {
    outFile->cd();
    hid->second->Write();
  }

  outFile->Write();
  delete outFile;
  outFile = 0;
}

void FillHist1D(const TString& histName, const Double_t& value, const double& weight) {
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight) {
  map<TString, TH2*>::iterator hid=m_Histos2D.find(histName);
  if (hid==m_Histos2D.end())
    cout << "%FillHist2D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}

void setWeight(const string& wFile) {

  ifstream file(wFile);
  string line;

  TString name = inName( inName.Last('/')+1, inName.Index('.')-inName.Last('/')-1 );
  if ( name.Contains("ttbar", TString::kIgnoreCase) ) {

    if      (topPtWeight != "NOMINAL") name += "_topPtWeight" + topPtWeight;
    else if (q2ttbar     != "NOMINAL") name += "_q2" + q2ttbar;
    else if (pdf         != "NOMINAL") name += "_pdf" + pdf;

    name.ReplaceAll("DOWN", "DN");
  }

  while (getline(file, line)) {

    if (line.length() > 0) {
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;

    TString dataset = line.substr(0, delim_pos).data();
    if ( name.EqualTo(dataset, TString::kIgnoreCase) ) {

      while (line.at(line.length()-1) == ' ') line.erase(line.length()-1, line.length());

      //weight is found in the last column
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      weight0 = stod(line);
      break;
    }
  }
  file.close();
}

void setPars(const string& parFile) {

  ifstream file(parFile);
  string line;

  while (getline(file, line)){

    if (line.length() > 0){
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;

    string var = line.substr(0, delim_pos);
    line.erase(0, delim_pos + 1);

    while (line.at(0) == ' ') line.erase(0, 1);
    while (line.at(line.length()-1) == ' ') line.erase(line.length()-1, line.length());

    if (var == "isMC"){
      if (line == "true") isMC = true;
      else isMC = false;
    }
    else if (var == "topPtWeight")   topPtWeight = line.data();
    else if (var == "jec")           jec = line.data();
    else if (var == "jer")           jer = line.data();
    else if (var == "btagSF")        btagSF = line.data();
    else if (var == "mistagSF")      mistagSF = line.data();
    else if (var == "pileup")        pileup = line.data();
    else if (var == "pdf")           pdf = line.data();
    else if (var == "q2ttbar")       q2ttbar = line.data();
    else if (var == "muTrigSys")     muTrigSys = line.data();
    else if (var == "muIdSys")       muIdSys = line.data();
    else if (var == "eleTrigSys")    eleTrigSys = line.data();
    else if (var == "eleIdSys")      eleIdSys = line.data();
    else if (var == "setSUMDRCut")   setSUMDRCut = line.data();
    else if (var == "inName")        inName = line.data();
    else if (var == "outName")       outName = line.data();
    else if (var == "muTrigSfName") {
      while ( (delim_pos = line.find(' ')) != -1) {
        int col = line.find(':');
        muTrigSfNames[line.substr(0, col).data()] = stod( line.substr(col+1, delim_pos-col-1) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      int col = line.find(':');
      muTrigSfNames[line.substr(0, col).data()] = stod( line.substr(col+1, delim_pos-col-1) );
    }
    else if (var == "muIdSfName") {
      while ( (delim_pos = line.find(' ')) != -1) {
        int col = line.find(':');
        muIdSfNames[line.substr(0, col).data()] = stod( line.substr(col+1, delim_pos-col-1) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      int col = line.find(':');
      muIdSfNames[line.substr(0, col).data()] = stod( line.substr(col+1, delim_pos-col-1) );
    }
    else if (var == "muTrackSfName") muTrackSfName = line.data();
    else if (var == "eTrigSfName")   eTrigSfName = line.data();
    else if (var == "eRecoSfName")   eRecoSfName = line.data();
    else if (var == "eIdSfName")     eIdSfName = line.data();
    else if (var == "btagName")      btagName = line.data();
    else if (var == "pileupName")    pileupName = line.data();
    else if (var == "channel")       channel = line;
    else if (var == "res_era")       res_era = line;
    else if (var == "eras") {
      while ( (delim_pos = line.find(' ')) != -1) {
        eras.push_back( line.substr(0, delim_pos) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      eras.push_back( line );
    }
    else if (var == "jet_type") jet_type = line;
  }
  file.close();
}

double rms_pm(const vector<float>& vec) {

  int size = vec.size();
  if (size == 0) return 0;
  double sum = 0;

  for (int i=0; i<size; i++) sum += ( vec[i]-1. )*( vec[i]-1. );

  return sqrt(sum / size);
}

bool newBTag( const float& coin, const float& pT, const int& flavor, const bool& oldBTag, TGraphAsymmErrors& g_eff, const TString& variation ) {
  double sf=0;

  //b or c jet
  //if ( abs(flavor) == 4 || abs(flavor) == 5 ) sf = 0.561694*((1.+(0.31439*pT))/(1.+(0.17756*pT)));   //medium SFs
  if ( abs(flavor) == 4 || abs(flavor) == 5 ) sf = 0.887973*((1.+(0.0523821*pT))/(1.+(0.0460876*pT))); //loose SFs, Run2016 BCDEFGH 
 
  //udsg
  //else sf = 1.06175-0.000462017*pT+1.02721e-06*pT*pT-4.95019e-10*pT*pT*pT;  //medium SFs
  //else sf = 1.15507+-0.00116691*pT+3.13873e-06*pT*pT+-2.14387e-09*pT*pT*pT; //loose SFs, Run2016 GH
  else sf = 1.12626+-0.000198068*pT+1.29872e-06*pT*pT+-1.00905e-09*pT*pT*pT;  //loose SFs, Run2016 BCDEF

  if (variation=="UP" || variation=="DOWN") {
    double sign = (variation=="UP") ? +1.: -1.;

    if (abs(flavor) == 4) {
       if (pT<30.)        sf += sign*0.063454590737819672 ;
       else if (pT<50.)   sf += sign*0.031410016119480133 ;
       else if (pT<70.)   sf += sign*0.02891194075345993  ;
       else if (pT<100.)  sf += sign*0.028121808543801308 ;
       else if (pT<140.)  sf += sign*0.027028990909457207 ;
       else if (pT<200.)  sf += sign*0.027206243947148323 ;
       else if (pT<300.)  sf += sign*0.033642303198575974 ;
       else if (pT<600.)  sf += sign*0.04273652657866478  ;
       else if (pT<1000.) sf += sign*0.054665762931108475 ;
       else               sf += sign*0.054665762931108475 * 2. ; // double syst. above 1 TeV
    }
    else if (abs(flavor) == 5) {
       if (pT<30.)        sf += sign*0.025381835177540779 ;
       else if (pT<50.)   sf += sign*0.012564006261527538 ;
       else if (pT<70.)   sf += sign*0.011564776301383972 ;
       else if (pT<100.)  sf += sign*0.011248723603785038 ;
       else if (pT<140.)  sf += sign*0.010811596177518368 ;
       else if (pT<200.)  sf += sign*0.010882497765123844 ;
       else if (pT<300.)  sf += sign*0.013456921093165874 ;
       else if (pT<600.)  sf += sign*0.017094610258936882 ;
       else if (pT<1000.) sf += sign*0.02186630479991436  ;
       else               sf += sign*0.02186630479991436 * 2. ; // double syst. above 1 TeV
    }
    else {
       if (pT<1000.) sf *= (1+sign*(0.100062-8.50875e-05*pT+4.8825e-08*pT*pT));
       else          sf *= (1+sign*(0.100062-8.50875e-05*pT+4.8825e-08*pT*pT) * 2.) ; // double syst. above 1 TeV
    }
  }

  if (sf == 1) return oldBTag; //no correction needed

  bool newBTag = oldBTag;

  if (sf > 1) {

    if( !oldBTag ) {

      float eff = g_eff.Eval( pT>1200?1200:pT );

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - sf) / (1.0 - (1.0/eff) );

      //upgrade to tagged
      if( coin < mistagPercent ) newBTag = true;
    }
  }
  else {
    //downgrade tagged to untagged
    if( oldBTag && coin > sf ) newBTag = false;
  }
  return newBTag;
}

