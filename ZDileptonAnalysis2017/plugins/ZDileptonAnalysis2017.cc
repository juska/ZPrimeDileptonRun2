// -*- C++ -*-
//
// Package:    analysis/ZDileptonAnalysis2017
// Class:      ZDileptonAnalysis2017
// 
//
//         Author:   Bahareh roozbahani Charles Harrington
//         Created:  Mon, 01 Oct 2018 16:26:57 GMT
//
//


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"   //For lumi block
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "parsePileUpJSON2.h"
#include <vector>
#include <algorithm>
#include <utility>
#include <TFile.h>
#include <TTree.h>
#include <TVector.h>
#include <TH1.h>
#include <TH2.h>


using namespace std;
using namespace edm;
bool sortLepPt(const pair<reco::CandidatePtr, char>& lep1, const pair<reco::CandidatePtr, char>& lep2){ return lep1.first->pt() > lep2.first->pt(); }

const int MAXJET = 50;
const int nFilters = 9;
const int MAXGEN = 20;
const int MAXLEP = 20;
const int METUNCERT = 4;
const int nTriggers = 10;

class ZDileptonAnalysis2017 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
      explicit ZDileptonAnalysis2017(const edm::ParameterSet&);
      ~ZDileptonAnalysis2017();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      double rms_pm(const vector<float>& vec) const;

      TFile* root_file;
      TTree* tree;
      TString fileName_;
      vector<int> totalEvts, dilep_cut, leppt_cut,dilepmass_cut;
      vector<int> jetpteta_cut, filter_failed, met_cut;
      vector<double> topPtWeightNOM, topPtWeightDN, q2UP, q2DN, q2UP_noTPW, q2DN_noTPW,  pdfUP, pdfDN, pdfUP_noTPW, pdfDN_noTPW;
      vector<float> wgt_env, wgt_rep;
      ULong64_t event;
      int run, lumi, bx;
      float rho, mu;
      int nPV, nPVall;
      char lep0flavor, lep1flavor;

      char jet_clean[MAXJET];
      int nJet, jet_hadflavor[MAXJET], jet_parflavor[MAXJET];
      float jet_pt[MAXJET], jet_eta[MAXJET], jet_phi[MAXJET], jet_mass[MAXJET], jet_area[MAXJET], jet_jec[MAXJET];
      float jet_btag[MAXJET], jet_btag_deepcsvprobb[MAXJET], jet_btag_deepcsvprobbb[MAXJET];
      float jet_nhf[MAXJET], jet_nef[MAXJET], jet_chf[MAXJET], jet_muf[MAXJET];
      float jet_elef[MAXJET], jet_numneutral[MAXJET], jet_chmult[MAXJET];

      int nGen;
      int gen_status[MAXGEN], gen_PID[MAXGEN], gen_mother0[MAXGEN], gen_mother1[MAXGEN], gen_index[MAXGEN];
      float gen_pt[MAXGEN], gen_mass[MAXGEN], gen_eta[MAXGEN], gen_phi[MAXGEN];
      int nGenJet;
      float genJet_pt[MAXJET], genJet_eta[MAXJET], genJet_phi[MAXJET], genJet_mass[MAXJET], genJet_area[MAXJET], genJet_nDaught[MAXJET];

      int nMuon;
      bool muon_isGlob[MAXLEP], muon_IsMediumID[MAXLEP], muon_IsTightID[MAXLEP];
      int muon_charge[MAXLEP];
      float muon_pt[MAXLEP], muon_eta[MAXLEP], muon_phi[MAXLEP], muon_D0[MAXLEP], muon_Dz[MAXLEP];
      float muon_chi2[MAXLEP], muon_tspm[MAXLEP], muon_kinkf[MAXLEP], muon_segcom[MAXLEP], muon_ftrackhits[MAXLEP];

      int nEle;
      int ele_charge[MAXLEP], ele_missinghits[MAXLEP];
      bool ele_LooseID[MAXLEP], ele_MediumID[MAXLEP], ele_TightID[MAXLEP]; //ele_passConv[MAXLEP];
      float ele_pt[MAXLEP], ele_eta[MAXLEP], ele_phi[MAXLEP], ele_D0[MAXLEP], ele_Dz[MAXLEP], ele_etaSupClust[MAXLEP];
      float ele_dPhiIn[MAXLEP], ele_sigmaIetaIeta[MAXLEP], ele_dEtaSeed[MAXLEP], ele_HE[MAXLEP], ele_rcpiwec[MAXLEP], ele_overEoverP[MAXLEP];

      float genmet_pt, genmet_px, genmet_py, genmet_sumet, genmet_phi;
      float met_pt, met_px, met_py, met_sumet, met_phi;
      int nMETUncert;
      float met_shiftedpx[METUNCERT], met_shiftedpy[METUNCERT];

      string metfilters[nFilters] = {
        "Flag_goodVertices",
        "Flag_globalSuperTightHalo2016Filter",
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter",
        "Flag_BadChargedCandidateFilter",
        "Flag_eeBadScFilter",
        "Flag_ecalBadCalibFilter"
      };
      vector<int> trig_prescale; vector<bool> trig_passed; vector<string> trig_name;

      string triggers[nTriggers] = {
        "HLT_Mu27_v",
        "HLT_Mu50_v",
        "HLT_Mu55_v",
        "HLT_TkMu100_v",
        "HLT_Mu37_TkMu27_v",
        "HLT_Mu27_Ele37_CaloIdL_MW_v",
        "HLT_Mu37_Ele27_CaloIdL_MW_v",
        "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele135_CaloIdVT_GsfTrkIdT_v",
        "HLT_DoubleEle33_CaloIdL_MW_v"
      };

      bool isMC_;
      edm::EDGetTokenT< edm::View<reco::GenParticle> > genParticleTag_;
      edm::EDGetTokenT<GenEventInfoProduct> genEventTag_;
      edm::EDGetTokenT<LHEEventProduct> extLHETag_;
      edm::EDGetTokenT< edm::View<PileupSummaryInfo> > muTag_;
      edm::EDGetTokenT< edm::View<pat::Muon> > muonTag_;
      edm::EDGetTokenT< edm::View<pat::Electron> > electronTag_;
      edm::EDGetTokenT< edm::View<pat::Jet> > jetTag_;
      double minLepPt_, minSubLepPt_, minDiLepMass_;
      string btag_, btag_deepcsvprobb_, btag_deepcsvprobbb_;
      double minLeadJetPt_;
      edm::EDGetTokenT<edm::TriggerResults> patTrgLabel_;
      edm::EDGetTokenT<double> rhoTag_;
      edm::EDGetTokenT< edm::View<reco::Vertex> > pvTag_;
      edm::EDGetTokenT< edm::View<reco::GenJet> > genJetTag_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleVetoIdMapToken_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleTightIdMapToken_;
      EffectiveAreas ele_areas_;
      edm::EDGetTokenT< edm::View<pat::MET> > metTag_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsTag_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> prescalesTag_;

      TH1F* fullMu = new TH1F("fullMu","fullMu",100,0,100);

      TH2F* ttbar_pt = new TH2F("ttbar_pt","ttbar_pt",250,0,5000,250,0,5000);
      TH2F* ttbar_pt2 = new TH2F("ttbar_pt2","ttbar_pt2",250,0,5000,250,0,5000);
      TH1F* deltat_pt = new TH1F("deltat_pt","deltat_pt",200,-1000,1000);
      TH1F* deltaTbar_pt = new TH1F("deltaTbar_pt","deltaTbar_pt",200,-1000,1000);

      TH1D* mass_totalEvts = new TH1D("mass_totalEvts","mass_totalEvts",500,0,10000);
      TH1D* mass_topPtWeightNOM = new TH1D("mass_topPtWeightNOM","mass_topPtWeightNOM",500,0,10000);
      TH1D* mass_topPtWeightDN = new TH1D("mass_topPtWeightDN","mass_topPtWeightDN",500,0,10000);
      TH1D* mass_q2UP = new TH1D("mass_q2UP","mass_q2UP",500,0,10000);
      TH1D* mass_q2DN = new TH1D("mass_q2DN","mass_q2DN",500,0,10000);
      TH1D* mass_q2UP_noTPW = new TH1D("mass_q2UP_noTPW","mass_q2UP_noTPW",500,0,10000);
      TH1D* mass_q2DN_noTPW = new TH1D("mass_q2DN_noTPW","mass_q2DN_noTPW",500,0,10000);
      TH1D* mass_pdfUP = new TH1D("mass_pdfUP","mass_pdfUP",500,0,10000);
      TH1D* mass_pdfDN = new TH1D("mass_pdfDN","mass_pdfDN",500,0,10000);
      TH1D* mass_pdfUP_noTPW = new TH1D("mass_pdfUP_noTPW","mass_pdfUP_noTPW",500,0,10000);
      TH1D* mass_pdfDN_noTPW = new TH1D("mass_pdfDN_noTPW","mass_pdfDN_noTPW",500,0,10000);
      TH2F* fullMu_mass = new TH2F("fullMu_mass","fullMu_mass",500,0,10000,100,0,100);
      TH1D* mass_dilep = new TH1D("mass_dilep","mass_dilep",500,0,10000);
      TH1D* mass_leppt = new TH1D("mass_leppt","mass_leppt",500,0,10000);
      TH1D* mass_dilepmass = new TH1D("mass_dilepmass","mass_dilepmass",500,0,10000);
      TH1D* mass_jetpteta = new TH1D("mass_jetpteta","mass_jetpteta",500,0,10000);
      TH1D* mass_met = new TH1D("mass_met","mass_met",500,0,10000);
};

ZDileptonAnalysis2017::ZDileptonAnalysis2017(const edm::ParameterSet& iConfig):
   ele_areas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )
{
  fileName_ = iConfig.getParameter<string>("fileName");
  isMC_ = iConfig.getParameter<bool>("isMC");
  genParticleTag_ = consumes< edm::View<reco::GenParticle> >( iConfig.getParameter<edm::InputTag>("genParticleTag") );
  genEventTag_ = consumes<GenEventInfoProduct>( iConfig.getParameter<edm::InputTag>("genEventTag") );
  extLHETag_ = consumes<LHEEventProduct>( iConfig.getParameter<edm::InputTag>("extLHETag") );
  muTag_ = consumes< edm::View<PileupSummaryInfo> >( iConfig.getParameter<edm::InputTag>("muTag") );
  muonTag_ = consumes< edm::View<pat::Muon> >( iConfig.getParameter<edm::InputTag>("muonTag") );
  electronTag_ = consumes< edm::View<pat::Electron> >( iConfig.getParameter<edm::InputTag>("electronTag") );
  minLepPt_ = iConfig.getParameter<double>("minLepPt");
  minSubLepPt_ = iConfig.getParameter<double>("minSubLepPt");
  minDiLepMass_ = iConfig.getParameter<double>("minDiLepMass");
  jetTag_ = consumes< edm::View<pat::Jet> >( iConfig.getParameter<edm::InputTag>("jetTag") );
  btag_ = iConfig.getParameter<string>("btag");
  btag_deepcsvprobb_ = iConfig.getParameter<string>("btag_deepcsvprobb");
  btag_deepcsvprobbb_ = iConfig.getParameter<string>("btag_deepcsvprobbb");
  minLeadJetPt_ = iConfig.getParameter<double>("minLeadJetPt");
  patTrgLabel_ = consumes<edm::TriggerResults>( iConfig.getParameter<edm::InputTag>("patTrgLabel") );
  rhoTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoTag") );
  pvTag_ = consumes< edm::View<reco::Vertex> >( iConfig.getParameter<edm::InputTag>("pvTag") );
  genJetTag_ = consumes< edm::View<reco::GenJet> >( iConfig.getParameter<edm::InputTag>("genJetTag") );
  eleVetoIdMapToken_ = consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("eleVetoIdMap") );
  eleLooseIdMapToken_ = consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("eleLooseIdMap") );
  eleMediumIdMapToken_ = consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("eleMediumIdMap") );
  eleTightIdMapToken_ = consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("eleTightIdMap") );
  metTag_ = consumes< edm::View<pat::MET> >( iConfig.getParameter<edm::InputTag>("metTag") );
  triggerResultsTag_ = consumes<edm::TriggerResults>( iConfig.getParameter<edm::InputTag>("triggerResultsTag") );
  prescalesTag_ = consumes<pat::PackedTriggerPrescales>( iConfig.getParameter<edm::InputTag>("prescalesTag") );
}

ZDileptonAnalysis2017::~ZDileptonAnalysis2017()
{
 

}

// ------------ method called once each job just before starting event loop  ------------
void 
ZDileptonAnalysis2017::beginJob()
{
  root_file = new TFile(fileName_, "RECREATE");
  tree = new TTree("T", "ZDilepton Tree");
  filter_failed.assign(nFilters+2, 0);
  totalEvts.assign(1, 0); dilep_cut.assign(1, 0); leppt_cut.assign(1, 0); jetpteta_cut.assign(1, 0); met_cut.assign(1, 0); dilepmass_cut.assign(1, 0);
  topPtWeightNOM.assign(1, 0.); topPtWeightDN.assign(1, 0.); pdfUP.assign(1, 0.); pdfDN.assign(1, 0.); q2UP.assign(1, 0.); q2DN.assign(1, 0.);
  pdfUP_noTPW.assign(1, 0.); pdfDN_noTPW.assign(1, 0.); q2UP_noTPW.assign(1, 0.); q2DN_noTPW.assign(1, 0.);

  tree->Branch("run", &run, "run/I");
  tree->Branch("lumi", &lumi, "lumi/I");
  tree->Branch("bx", &bx, "bx/I");
  tree->Branch("event", &event, "event/l");
  tree->Branch("rho", &rho, "rho/F");
  tree->Branch("mu", &mu, "mu/F");
  tree->Branch("nPV", &nPV, "nPV/I");
  tree->Branch("nPVall", &nPVall, "nPVall/I");
  tree->Branch("lep0flavor", &lep0flavor, "lep0flavor/B");
  tree->Branch("lep1flavor", &lep1flavor, "lep1flavor/B");
  tree->Branch("nJet", &nJet, "nJet/I");
  tree->Branch("jet_clean", jet_clean, "jet_clean[nJet]/B");
  tree->Branch("jet_pt", jet_pt, "jet_pt[nJet]/F");
  tree->Branch("jet_eta", jet_eta, "jet_eta[nJet]/F");
  tree->Branch("jet_phi", jet_phi, "jet_phi[nJet]/F");
  tree->Branch("jet_mass", jet_mass, "jet_mass[nJet]/F");
  tree->Branch("jet_area", jet_area, "jet_area[nJet]/F");
  tree->Branch("jet_jec", jet_jec, "jet_jec[nJet]/F");
  tree->Branch("jet_btag", jet_btag, "jet_btag[nJet]/F");
  tree->Branch("jet_btag_deepcsvprobb", jet_btag_deepcsvprobb, "jet_btag_deepcsvprobb[nJet]/F");
  tree->Branch("jet_btag_deepcsvprobbb", jet_btag_deepcsvprobbb, "jet_btag_deepcsvprobbb[nJet]/F");
  tree->Branch("nMuon", &nMuon, "nMuon/I");
  tree->Branch("muon_charge", muon_charge, "muon_charge[nMuon]/I");
  tree->Branch("muon_isGlob", muon_isGlob, "muon_isGlob[nMuon]/O");
  tree->Branch("muon_pt", muon_pt, "muon_pt[nMuon]/F");
  tree->Branch("muon_eta", muon_eta, "muon_eta[nMuon]/F");
  tree->Branch("muon_phi", muon_phi, "muon_phi[nMuon]/F");
  tree->Branch("muon_D0", muon_D0, "muon_D0[nMuon]/F");
  tree->Branch("muon_Dz", muon_Dz, "muon_Dz[nMuon]/F");
  tree->Branch("muon_chi2", muon_chi2, "muon_chi2[nMuon]/F");
  tree->Branch("muon_tspm", muon_tspm, "muon_tspm[nMuon]/F");
  tree->Branch("muon_kinkf", muon_kinkf, "muon_kinkf[nMuon]/F");
  tree->Branch("muon_ftrackhits", muon_ftrackhits, "muon_ftrackhits[nMuon]/F");
  tree->Branch("muon_segcom", muon_segcom, "muon_segcom[nMuon]/F");
  tree->Branch("muon_IsMediumID", muon_IsMediumID, "muon_IsMediumID[nMuon]/O");
  tree->Branch("muon_IsTightID", muon_IsTightID, "muon_IsTightID[nMuon]/O");
  tree->Branch("nEle", &nEle, "nEle/I");
  tree->Branch("ele_charge", ele_charge, "ele_charge[nEle]/I");
  tree->Branch("ele_pt", ele_pt, "ele_pt[nEle]/F");
  tree->Branch("ele_eta", ele_eta, "ele_eta[nEle]/F");
  tree->Branch("ele_phi", ele_phi, "ele_phi[nEle]/F");
  tree->Branch("ele_LooseID", ele_LooseID , "ele_LooseID[nEle]/O");
  tree->Branch("ele_MediumID", ele_MediumID , "ele_MediumID[nEle]/O");
  tree->Branch("ele_TightID", ele_TightID , "ele_TightID[nEle]/O");
  tree->Branch("ele_D0", ele_D0, "ele_D0[nEle]/F");
  tree->Branch("ele_Dz", ele_Dz, "ele_Dz[nEle]/F");
  tree->Branch("ele_sigmaIetaIeta", ele_sigmaIetaIeta, "ele_sigmaIetaIeta[nEle]/F");
  tree->Branch("ele_dEtaSeed", ele_dEtaSeed, "ele_dEtaSeed[nEle]/F");
  tree->Branch("ele_dPhiIn", ele_dPhiIn, "ele_dPhiIn[nEle]/F");
  tree->Branch("ele_etaSupClust", ele_etaSupClust, "ele_etaSupClust[nEle]/F");
  tree->Branch("ele_overEoverP", ele_overEoverP, "ele_overEoverP[nEle]/F");
  tree->Branch("ele_HE", ele_HE, "ele_HE[nEle]/F");
  tree->Branch("ele_rcpiwec", ele_rcpiwec, "ele_rcpiwec[nEle]/F");
  tree->Branch("ele_missinghits", ele_missinghits, "ele_missinghits[nEle]/I");
  tree->Branch("met_pt", &met_pt, "met_pt/F");
  tree->Branch("met_px", &met_px, "met_px/F");
  tree->Branch("met_py", &met_py, "met_py/F");
  tree->Branch("met_sumet", &met_sumet, "met_sumet/F");
  tree->Branch("met_phi", &met_phi, "met_phi/F");
  tree->Branch("nMETUncert", &nMETUncert, "nMETUncert/I");
  tree->Branch("met_shiftedpx", met_shiftedpx, "met_shiftedpx[nMETUncert]/F");
  tree->Branch("met_shiftedpy", met_shiftedpy, "met_shiftedpy[nMETUncert]/F");
  tree->Branch("trig_prescale", "std::vector<int>", &trig_prescale);
  tree->Branch("trig_passed", "std::vector<bool>", &trig_passed);
  tree->Branch("trig_name", "std::vector<string>", &trig_name);
  if(isMC_) {
    tree->Branch("wgt_env", "std::vector<float>", &wgt_env);
    tree->Branch("wgt_rep", "std::vector<float>", &wgt_rep);
    tree->Branch("jet_hadflavor", jet_hadflavor, "jet_hadflavor[nJet]/I");
    tree->Branch("jet_parflavor", jet_parflavor, "jet_parflavor[nJet]/I");

    tree->Branch("nGen", &nGen, "nGen/I");
    tree->Branch("gen_status",  gen_status, "gen_status[nGen]/I");
    tree->Branch("gen_PID",  gen_PID, "gen_PID[nGen]/I");
    tree->Branch("gen_pt", gen_pt, "gen_pt[nGen]/F");
    tree->Branch("gen_mass", gen_mass, "gen_mass[nGen]/F");
    tree->Branch("gen_eta", gen_eta, "gen_eta[nGen]/F");
    tree->Branch("gen_phi", gen_phi, "gen_phi[nGen]/F");
    tree->Branch("gen_index",  gen_index, "gen_index[nGen]/I");
    tree->Branch("gen_mother0", gen_mother0, "gen_mother0[nGen]/I");
    tree->Branch("gen_mother1", gen_mother1, "gen_mother1[nGen]/I");

    tree->Branch("nGenJet", &nGenJet, "nGenJet/I");
    tree->Branch("genJet_pt", genJet_pt, "genJet_pt[nGenJet]/F");
    tree->Branch("genJet_eta", genJet_eta, "genJet_eta[nGenJet]/F");
    tree->Branch("genJet_phi", genJet_phi, "genJet_phi[nGenJet]/F");
    tree->Branch("genJet_mass", genJet_mass, "genJet_mass[nGenJet]/F");
    tree->Branch("genJet_area",  genJet_area, "genJet_area[nGenJet]/F");
    tree->Branch("genJet_nDaught", genJet_nDaught, "genJet_nDaught[nGenJet]/F");
    tree->Branch("genmet_pt", &genmet_pt, "genmet_pt/F");
    tree->Branch("genmet_px", &genmet_px, "genmet_px/F");
    tree->Branch("genmet_py", &genmet_py, "genmet_py/F");
    tree->Branch("genmet_sumet", &genmet_sumet, "genmet_sumet/F");
    tree->Branch("genmet_phi", &genmet_phi, "genmet_phi/F");
  }

  else{
    parsePileUpJSON2();
  }
}

// ------------ method called for each event  ------------
void
ZDileptonAnalysis2017::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  totalEvts[0]++;

  //------------ Event Info ------------//

  run = int(iEvent.id().run());
  lumi = int(iEvent.getLuminosityBlock().luminosityBlock());
  bx = iEvent.bunchCrossing();
  event = iEvent.id().event();
  cout << " event " << event << endl;
  double mass_ttbar = 0;

 //------------ MC :topPt, pdf, and q2 ------------//
  if (isMC_) {
    edm::Handle< edm::View<reco::GenParticle> > genParticles;
    iEvent.getByToken(genParticleTag_, genParticles);

    reco::Candidate::LorentzVector t_p4, tbar_p4;
    float t_pt2=-1, tbar_pt2=-1;
    for (int i=0, n=genParticles->size(); i<n; i++) {
      const reco::GenParticle& p = genParticles->at(i);

      int id = p.pdgId();
      int status = p.status();
      int nDaught = p.numberOfDaughters();

      if (id==6 && 20<=status && status<30) { t_p4 = p.p4(); cout << i << "\t t first " << t_p4.Pt() << endl;}
      else if (id==-6 && 20<=status && status<30) tbar_p4 = p.p4();

      else if (fabs(id)==6 && nDaught==2) {   //last t's
        const reco::GenParticle& daught0 = genParticles->at( p.daughterRef(0).key() );

        if ( fabs(daught0.pdgId())==5 || fabs(daught0.pdgId())==24 ) {
          if (id==6) { t_pt2 = p.pt(); cout << i << "\t t last " << t_pt2 << endl;}
          else       tbar_pt2 = p.pt();
        }
      }
    }
    float t_pt = t_p4.Pt(), tbar_pt = tbar_p4.Pt();
    mass_ttbar = (t_p4+tbar_p4).M() > 10000 ? 9999 : (t_p4+tbar_p4).M();

    ttbar_pt->Fill(t_pt, tbar_pt);
    ttbar_pt2->Fill(t_pt2, tbar_pt2);
    deltat_pt->Fill(t_pt-t_pt2);
    deltaTbar_pt->Fill(tbar_pt-tbar_pt2);
    mass_totalEvts->Fill(mass_ttbar);

    double wgt_topPtWeightNOM = sqrt( exp(0.0615-0.0005*t_pt2) * exp(0.0615-0.0005*tbar_pt2) );
    double wgt_topPtWeightDN = exp(0.0615-0.0005*t_pt2) * exp(0.0615-0.0005*tbar_pt2);

    topPtWeightNOM[0] += wgt_topPtWeightNOM;
    topPtWeightDN[0] += wgt_topPtWeightDN;
    mass_topPtWeightNOM->Fill( mass_ttbar, wgt_topPtWeightNOM );
    mass_topPtWeightDN->Fill( mass_ttbar, wgt_topPtWeightDN );

    wgt_env.clear(); wgt_rep.clear();
    edm::Handle<GenEventInfoProduct> genEventHandle;
    //nominal weight = 1, twiki says this must be activated
    if ( iEvent.getByToken(genEventTag_, genEventHandle) ) genEventHandle->weight(); // why we dont use this?
    edm::Handle<LHEEventProduct> lheEvtProduct;
    if ( iEvent.getByToken(extLHETag_ , lheEvtProduct) ) {
      float wgt_denom = lheEvtProduct->weights()[0].wgt;
      for (int i=0; i<9; i++) { // push back the seven envelope weights
        if (i == 5 || i == 7) continue;
        wgt_env.push_back( lheEvtProduct->weights()[i].wgt/wgt_denom );
      }
      double wgt_q2UP = TMath::MaxElement(wgt_env.size(), &wgt_env[0]);
      double wgt_q2DN = TMath::MinElement(wgt_env.size(), &wgt_env[0]);
      q2UP[0] += wgt_q2UP * wgt_topPtWeightNOM;  q2UP_noTPW[0] += wgt_q2UP;
      q2DN[0] += wgt_q2DN * wgt_topPtWeightNOM; q2DN_noTPW[0] += wgt_q2DN;
      mass_q2UP->Fill( mass_ttbar, wgt_q2UP * wgt_topPtWeightNOM );  mass_q2UP_noTPW->Fill( mass_ttbar, wgt_q2UP );
      mass_q2DN->Fill( mass_ttbar, wgt_q2DN * wgt_topPtWeightNOM ); mass_q2DN_noTPW->Fill( mass_ttbar, wgt_q2DN );


      // push back 100 weight replicas
      for (int i=9, n=9+100; i<n; i++) wgt_rep.push_back( lheEvtProduct->weights()[i].wgt/wgt_denom );

      vector<float> pdf_plus, pdf_minus;
      for (unsigned int i=0, n=wgt_rep.size(); i<n; i++) {
        if (wgt_rep[i] >= 1.) pdf_plus.push_back(wgt_rep[i]);
        else                  pdf_minus.push_back(wgt_rep[i]);
      }
      double wgt_pdfUP = 1. + rms_pm(pdf_plus);
      double wgt_pdfDN = 1. - rms_pm(pdf_minus);
     pdfUP[0] += wgt_pdfUP * wgt_topPtWeightNOM;  pdfUP_noTPW[0] += wgt_pdfUP;
     pdfDN[0] += wgt_pdfDN * wgt_topPtWeightNOM; pdfDN_noTPW[0] += wgt_pdfDN;

     mass_pdfUP->Fill( mass_ttbar, wgt_pdfUP * wgt_topPtWeightNOM );  mass_pdfUP_noTPW->Fill( mass_ttbar, wgt_pdfUP );
     mass_pdfDN->Fill( mass_ttbar, wgt_pdfDN * wgt_topPtWeightNOM ); mass_pdfDN_noTPW->Fill( mass_ttbar, wgt_pdfDN );
    }


    edm::Handle< edm::View<PileupSummaryInfo> > pileups;
    iEvent.getByToken(muTag_, pileups);
    mu = pileups->at(1).getTrueNumInteractions();
    fullMu_mass->Fill(mass_ttbar, mu);
  }
  else {
    mu = getAvgPU( run, lumi );
  }
  fullMu->Fill(mu);

  //------------ Lepton Pt Filter ------------//
  edm::Handle< edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonTag_, muons);
  edm::Handle< edm::View<pat::Electron> > electrons;
  iEvent.getByToken(electronTag_, electrons);
  edm::Handle<edm::ValueMap<bool> > Veto_id_decisions;
  iEvent.getByToken(eleVetoIdMapToken_ ,Veto_id_decisions);

  vector<pair<reco::CandidatePtr, char> > leps;

  for (int i=0, n=muons->size(); i<n; i++) {
    if (!muons->at(i).isLooseMuon()) continue;
    leps.push_back( make_pair(muons->ptrAt(i),'m') );
  }
  for (int i=0, n=electrons->size(); i<n; i++) {
    const Ptr<pat::Electron> elPtr(electrons, i);
    if (!(*Veto_id_decisions)[elPtr]) continue;
    leps.push_back( make_pair(electrons->ptrAt(i),'e') );
  }


  if (leps.size() < 2) return;
  dilep_cut[0]++;
  mass_dilep->Fill(mass_ttbar);

  sort(leps.begin(), leps.end(), sortLepPt);
  if (leps[0].first->pt() < minLepPt_ || leps[1].first->pt() < minSubLepPt_) return;
  leppt_cut[0]++;
  mass_leppt->Fill(mass_ttbar);

  lep0flavor = leps[0].second;
  lep1flavor = leps[1].second;

  if (lep0flavor == lep1flavor) {
    bool dilepmass_flag = false;
    for (int i=0, n=leps.size(); i<n && !dilepmass_flag; i++) {
      char flavi = leps[i].second;
      reco::Candidate::LorentzVector vlepi = leps[i].first->p4();

      for (int j=i+1; j<n && !dilepmass_flag; j++) {
        char flavj = leps[j].second;

        if (flavi==flavj) {
          reco::Candidate::LorentzVector vlepj = leps[j].first->p4();
          if ( (vlepi+vlepj).M() > minDiLepMass_ ) dilepmass_flag = true;
        }
      }
    }
    if (!dilepmass_flag) return;
  }
  dilepmass_cut[0]++;
  mass_dilepmass->Fill(mass_ttbar);

  vector<reco::CandidatePtr> lep0Sources, lep1Sources;
  for (unsigned int i=0, n=leps[0].first->numberOfSourceCandidatePtrs(); i<n; i++) lep0Sources.push_back(leps[0].first->sourceCandidatePtr(i));
  for (unsigned int i=0, n=leps[1].first->numberOfSourceCandidatePtrs(); i<n; i++) lep1Sources.push_back(leps[1].first->sourceCandidatePtr(i));


  //------------ Jets ------------//

  edm::Handle< edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetTag_, jets);

  int jetsize = jets->size();
  nJet = jetsize>MAXJET ? MAXJET : jetsize;
  int leadJetPt_flag = false;
  for (int i=0; i<nJet; i++){
    const pat::Jet& jet = jets->at(i).correctedJet(0);
    reco::Candidate::LorentzVector jet_p4 = jet.p4();
    jet_clean[i] = 'n';
    if ( reco::deltaR(jet, *leps[0].first) < 0.4 || reco::deltaR(jet, *leps[1].first) < 0.4 ){
      const vector<reco::CandidatePtr> & dvec = jet.daughterPtrVector();
      for (vector<reco::CandidatePtr>::const_iterator i_d = dvec.begin(); i_d != dvec.end(); ++i_d){
        if ( find(lep0Sources.begin(), lep0Sources.end(), *i_d ) != lep0Sources.end() ) {
          jet_p4 -= (*i_d)->p4();
          if (jet_clean[i] == 'n') jet_clean[i] = 'l';
          else if (jet_clean[i] == 's') jet_clean[i] = 'b';
        } 
        else if ( find(lep1Sources.begin(), lep1Sources.end(), *i_d ) != lep1Sources.end() ) {
          jet_p4 -= (*i_d)->p4();
          if (jet_clean[i] == 'n') jet_clean[i] = 's';
          else if (jet_clean[i] == 'l') jet_clean[i] = 'b';
        } 

      }
    }

    jet_pt[i] = jet_p4.Pt();
    jet_eta[i] = jet_p4.Eta();
    jet_phi[i] = jet_p4.Phi();
    jet_mass[i] = jet_p4.M();
    jet_area[i] = jet.jetArea();
    jet_jec[i] = jets->at(i).jecFactor(0);

    jet_nhf[i] = jet.neutralHadronEnergyFraction();
    jet_nef[i] = jet.neutralEmEnergyFraction();
    jet_chf[i] = jet.chargedHadronEnergyFraction();
    jet_muf[i] = jet.muonEnergyFraction();
    jet_elef[i] = jet.chargedEmEnergyFraction();
    jet_numneutral[i] =jet.neutralMultiplicity();
    jet_chmult[i] = jet.chargedMultiplicity();
    jet_btag[i] = jet.bDiscriminator(btag_);
    jet_btag_deepcsvprobb[i] = jet.bDiscriminator(btag_deepcsvprobb_);
    jet_btag_deepcsvprobbb[i] = jet.bDiscriminator(btag_deepcsvprobbb_);


    if (isMC_) {
      jet_hadflavor[i] = jet.hadronFlavour();
      jet_parflavor[i] = jet.partonFlavour();
    }

    if ( jet_pt[i]>minLeadJetPt_ && fabs(jet_eta[i])<2.5 ) leadJetPt_flag = true;

  }
  if (!leadJetPt_flag) return;
  jetpteta_cut[0]++;
  mass_jetpteta->Fill(mass_ttbar);


 //------------ MET Filters ------------//

  Handle<edm::TriggerResults> patFilterHandle;
  iEvent.getByToken(patTrgLabel_, patFilterHandle);
  if (patFilterHandle.isValid()){
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*patFilterHandle);
    for (int i=0; i<nFilters; i++){
      if(  isMC_ && metfilters[i] == "Flag_eeBadScFilter" ) continue;
      const unsigned int trig = trigNames.triggerIndex(metfilters[i]);
      if (trig != trigNames.size()){
        if (!patFilterHandle->accept(trig)){
          filter_failed[i]++;
          return;
        }
      }
      else cout << metfilters[i] << " not found." << endl;
    }

  }

  met_cut[0]++;
  mass_met->Fill(mass_ttbar);


  //------------ Rho ------------//

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoTag_, rhoHandle);
  rho = *rhoHandle;


  //------------ Primary Vertices ------------//

  edm::Handle< edm::View<reco::Vertex> > primaryVertices;
  iEvent.getByToken(pvTag_, primaryVertices);

  reco::Vertex pvtx = primaryVertices->at(0);

  nPVall = primaryVertices->size();
  nPV = 0;

  for (int i=0; i<nPVall; i++) {
    const reco::Vertex& pv = primaryVertices->at(i);
    if( !pv.isFake() && pv.ndof() > 4 && abs(pv.z()) <= 24 && pv.position().rho() <= 2 )
      nPV++;
  }

//-------------- Generated Particles -------------//

  if(isMC_) {
    edm::Handle< edm::View<reco::GenParticle> > genParticles;
    iEvent.getByToken(genParticleTag_, genParticles);
    vector<pair<reco::GenParticle, int> > reducedGens;
    for (int i=0, n=genParticles->size(); i<n; i++) {
      const reco::GenParticle& p = genParticles->at(i);
      int id = p.pdgId();
      int status = p.status();
      int nDaught = p.numberOfDaughters();
      if (id>1000000)  //Z'
        reducedGens.push_back(make_pair(p,i));
      if (fabs(id)==6 && 20<=status && status<30)  //first t's
        reducedGens.push_back(make_pair(p,i));

      else if (fabs(id)==6 && nDaught==2) {   //last t's
        const reco::GenParticle& daught0 = genParticles->at( p.daughterRef(0).key() );
        if ( fabs(daught0.pdgId())==5 || fabs(daught0.pdgId())==24 ) {
          reducedGens.push_back(make_pair(p,i));
          reducedGens.push_back(make_pair( daught0,p.daughterRef(0).key()) ); //b or W (first)
          reducedGens.push_back(make_pair( genParticles->at(p.daughterRef(1).key()),p.daughterRef(1).key()) ); //b or W (first)
        }
      }


      else if (fabs(id)==24 && nDaught==2) {   //last W's
        reducedGens.push_back(make_pair(p,i));
        reducedGens.push_back(make_pair( genParticles->at(p.daughterRef(0).key()),p.daughterRef(0).key()) ); //q or lep
        reducedGens.push_back(make_pair( genParticles->at(p.daughterRef(1).key()),p.daughterRef(1).key()) ); //q or lep
      }
    }
    nGen = reducedGens.size();
    for (int i=0; i<nGen; i++) {
      const reco::GenParticle& p = reducedGens[i].first;
      gen_status[i] = p.status();
      gen_PID[i] = p.pdgId();
      gen_pt[i]  = p.pt();
      gen_mass[i] = p.mass();
      gen_eta[i] = p.eta();
      gen_phi[i] = p.phi();
      gen_index[i] = reducedGens[i].second;
      if (p.numberOfMothers() > 0) gen_mother0[i] = p.motherRef(0).key();
      else gen_mother0[i] = -1;
      if (p.numberOfMothers() > 1) gen_mother1[i] = p.motherRef(1).key();
      else gen_mother1[i] = -1;

    }

    edm::Handle< edm::View<reco::GenJet> > genJets;
    iEvent.getByToken(genJetTag_, genJets);

    nGenJet = 0;

    for (int i=0, n=genJets->size(); i<n; i++){
      if (genJets->at(i).pt() < 15) continue;
      const reco::GenJet& jet = genJets->at(i);
      genJet_pt[nGenJet]  = jet.pt();
      genJet_eta[nGenJet] = jet.eta();
      genJet_phi[nGenJet] = jet.phi();
      genJet_mass[nGenJet] = jet.mass();
      genJet_area[nGenJet] = jet.jetArea();
      genJet_nDaught[nGenJet] = jet.numberOfDaughters();

      nGenJet++;
    }
  }

  //--------------Muons-------------//
  nMuon = 0;
  for (int i=0, n=muons->size(); i<n; i++){
    //save only loose muons above 15 gev
    if ( muons->at(i).pt() < 15 || !muons->at(i).isLooseMuon() ) continue;
    const pat::Muon& muon = muons->at(i);
    muon_isGlob[nMuon] = muon.isGlobalMuon();
    cout << muon_isGlob[nMuon] << endl;
    muon_charge[nMuon] = muon.charge();
    muon_pt[nMuon] = muon.pt();
    muon_eta[nMuon] = muon.eta();
    muon_phi[nMuon] = muon.phi();

    muon_D0[nMuon] = muon.muonBestTrack()->dxy(pvtx.position());
    muon_Dz[nMuon] = muon.muonBestTrack()->dz(pvtx.position());
    muon_tspm[nMuon] = muon.combinedQuality().chi2LocalPosition;
    muon_kinkf[nMuon] = muon.combinedQuality().trkKink;
    muon_ftrackhits[nMuon] = muon.innerTrack()->validFraction();

    muon_segcom[nMuon] = muon::segmentCompatibility(muon);

    if (muon_isGlob[nMuon]) muon_chi2[nMuon] = muon.globalTrack()->normalizedChi2();
    else                    muon_chi2[nMuon] = 99;

    muon_IsMediumID[nMuon] = muon.isMediumMuon();
    muon_IsTightID[nMuon] = muon.isTightMuon(pvtx);   

    nMuon++;
  }


  //------------ Electrons ------------//

  edm::Handle<edm::ValueMap<bool> > Loose_id_decisions;
  iEvent.getByToken(eleLooseIdMapToken_, Loose_id_decisions);
  edm::Handle<edm::ValueMap<bool> > Medium_id_decisions;
  iEvent.getByToken(eleMediumIdMapToken_, Medium_id_decisions);
  edm::Handle<edm::ValueMap<bool> > Tight_id_decisions;
  iEvent.getByToken(eleTightIdMapToken_, Tight_id_decisions);
  nEle = 0;
  for (int i=0, n=electrons->size(); i<n; i++){
    const pat::Electron& ele = electrons->at(i);
    const Ptr<pat::Electron> elPtr(electrons, i);
   //save only veto electrons (without isocut) above 15 gev
    if (ele.pt() < 15 || !(*Veto_id_decisions)[elPtr]) continue;
    ele_charge[nEle] = ele.charge();
    ele_pt[nEle] = ele.pt();
    ele_eta[nEle] = ele.eta();
    ele_phi[nEle] = ele.phi();
    ele_D0[nEle] = ele.gsfTrack()->dxy(pvtx.position());
    ele_Dz[nEle] = ele.gsfTrack()->dz(pvtx.position());
    ele_sigmaIetaIeta[nEle] = ele.full5x5_sigmaIetaIeta();

    if (ele.superCluster().isNonnull()) ele_etaSupClust[nEle] = ele.superCluster()->eta();
    else ele_etaSupClust[nEle] = -99;
    if ( ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull() )
      ele_dEtaSeed[nEle] = ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta();
    else ele_dEtaSeed[nEle] = -99;

    ele_dPhiIn[nEle] = ele.deltaPhiSuperClusterTrackAtVtx();
    ele_overEoverP[nEle] = fabs(1.0 - ele.eSuperClusterOverP()) / ele.ecalEnergy();

    ele_HE[nEle] = ele.hadronicOverEm();
    constexpr reco::HitPattern::HitCategory missingHitType = reco::HitPattern::MISSING_INNER_HITS;
    ele_missinghits[nEle] = ele.gsfTrack()->hitPattern().numberOfAllHits(missingHitType); // This  has changed to numberOfAllHits //check

    reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
    float eA = ele_areas_.getEffectiveArea( fabs(ele_etaSupClust[nEle]) );
    ele_rcpiwec[nEle] = ( pfIso.sumChargedHadronPt + max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho) ) / ele_pt[nEle];
    ele_LooseID[nEle]  = (*Loose_id_decisions)[elPtr];
    ele_MediumID[nEle] = (*Medium_id_decisions)[elPtr];
    ele_TightID[nEle]  = (*Tight_id_decisions)[elPtr];
    nEle++;
  }

  //------------ MET ------------//

  edm::Handle< edm::View<pat::MET> > mets;
  iEvent.getByToken(metTag_, mets);
  const pat::MET& met = mets->at(0);
  if (isMC_){
    const reco::GenMET *genmet = met.genMET();
    genmet_pt = genmet->pt();
    genmet_px = genmet->px();
    genmet_py = genmet->py();
    genmet_sumet = genmet->sumEt();
    genmet_phi = genmet->phi();
  }

  met_pt = met.uncorPt();
  met_px = met.uncorPx();
  met_py = met.uncorPy();
  met_sumet = met.uncorSumEt();
  met_phi = met.uncorPhi();

  nMETUncert = METUNCERT;

  for (int i=0; i<nMETUncert; i++){
    pat::MET::METUncertainty uncert = static_cast<pat::MET::METUncertainty>(i);
    pat::MET::METCorrectionLevel level = pat::MET::Type1;
    met_shiftedpx[i] = met.shiftedPx(uncert, level);
    met_shiftedpy[i] = met.shiftedPy(uncert, level);
  }

  //------------ Triggers ------------//

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsTag_, triggerResults);
  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(prescalesTag_, triggerPrescales);

  if (triggerResults.isValid()) {
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);   
    trig_prescale.assign(nTriggers, 0); trig_passed.assign(nTriggers, false); trig_name.assign(nTriggers, "");
    for (int i=0; i<nTriggers; i++) {
      string myTrigger = triggers[i], name;
      int len = myTrigger.length(), index = -1;

      for (int i=0, n=triggerResults->size(); i<n; i++){
        name = triggerNames.triggerName(i);
        if (name.substr(0, len) == myTrigger) {index = i; cout << "run \t" << run << "\t event \t" << "\texisting trigger\t" << name << endl; break;}
      }
      if (index == -1) continue;

      trig_prescale[i] = triggerPrescales->getPrescaleForIndex(index);
      trig_passed[i] = triggerResults->accept(index);
      trig_name[i] = name;
    }

  }

  //------------ Fill Tree ------------//
  tree->Fill();

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZDileptonAnalysis2017::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZDileptonAnalysis2017::endJob() 
{
  if (root_file !=0) {
    root_file->cd();
    fullMu->Write();
    if (isMC_) {
      fullMu_mass->Write();
      mass_totalEvts->Write();
      mass_dilep->Write();
      mass_leppt->Write();
      mass_dilepmass->Write();
      mass_jetpteta->Write();
      mass_met->Write();

      mass_topPtWeightNOM->Write();
      mass_topPtWeightDN->Write();
      mass_q2UP->Write();
      mass_q2DN->Write();
      mass_q2UP_noTPW->Write();
      mass_q2DN_noTPW->Write();
      mass_pdfUP->Write();
      mass_pdfDN->Write();
      mass_pdfUP_noTPW->Write();
      mass_pdfDN_noTPW->Write();

      ttbar_pt->Write();
      ttbar_pt2->Write();
      deltat_pt->Write();
      deltaTbar_pt->Write();

      root_file->WriteObject(&topPtWeightNOM, "topPtWeightNOM");
      root_file->WriteObject(&topPtWeightDN, "topPtWeightDN");
      root_file->WriteObject(&q2UP, "q2UP");
      root_file->WriteObject(&q2DN, "q2DN");
      root_file->WriteObject(&q2UP_noTPW, "q2UP_noTPW");
      root_file->WriteObject(&q2DN_noTPW, "q2DN_noTPW");
      root_file->WriteObject(&pdfUP, "pdfUP");
      root_file->WriteObject(&pdfDN, "pdfDN");
      root_file->WriteObject(&pdfUP_noTPW, "pdfUP_noTPW");
      root_file->WriteObject(&pdfDN_noTPW, "pdfDN_noTPW");
    }
    root_file->WriteObject(&totalEvts, "totalEvts");
    root_file->WriteObject(&dilep_cut, "dilep_cut");
    root_file->WriteObject(&leppt_cut, "leppt_cut");
    root_file->WriteObject(&dilepmass_cut, "dilepmass_cut");
    root_file->WriteObject(&jetpteta_cut, "jetpteta_cut");
    root_file->WriteObject(&filter_failed, "filter_failed");
    root_file->WriteObject(&met_cut, "met_cut");


    root_file->Write();
    delete root_file;
    root_file = 0;
  }
}

double ZDileptonAnalysis2017::rms_pm(const vector<float>& vec) const {

  int size = vec.size();
  if (size == 0) return 0;
  double sum = 0;

  for (int i=0; i<size; i++) sum += ( vec[i]-1. )*( vec[i]-1. );

  return sqrt(sum / size);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZDileptonAnalysis2017);
