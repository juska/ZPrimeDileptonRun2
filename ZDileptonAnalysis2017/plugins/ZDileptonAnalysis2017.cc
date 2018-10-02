// -*- C++ -*-
//
// Package:    analysis/ZDileptonAnalysis2017
// Class:      ZDileptonAnalysis2017
// 
//
// Original Author:  bahareh roozbahani
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


class ZDileptonAnalysis2017 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
      explicit ZDileptonAnalysis2017(const edm::ParameterSet&);
      ~ZDileptonAnalysis2017();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

    TFile* root_file;
    TTree* tree;
    TString fileName_;
    vector<int> totalEvts;
    ULong64_t event;
    int run, lumi, bx;
    float rho, mu;
    int nPV, nPVall;
    bool isMC_;
    edm::EDGetTokenT< edm::View<reco::GenParticle> > genParticleTag_;
    edm::EDGetTokenT< edm::View<PileupSummaryInfo> > muTag_;
    edm::EDGetTokenT< edm::View<pat::Muon> > muonTag_;
    edm::EDGetTokenT< edm::View<pat::Electron> > electronTag_;
    edm::EDGetTokenT< edm::ValueMap<bool> > eleVetoIdMapToken_;


    TH1D* mass_totalEvts = new TH1D("mass_totalEvts","mass_totalEvts",500,0,5000);
    TH2F* ttbar_pt = new TH2F("ttbar_pt","ttbar_pt",50,0,1000,50,0,1000);
    TH2F* ttbar_pt2 = new TH2F("ttbar_pt2","ttbar_pt2",50,0,1000,50,0,1000);
    TH1F* deltat_pt = new TH1F("deltat_pt","deltat_pt",100,-500,500);
    TH1F* deltaTbar_pt = new TH1F("deltaTbar_pt","deltaTbar_pt",100,-500,500);
    TH1F* fullMu = new TH1F("fullMu","fullMu",100,0,100);
    TH2F* fullMu_mass = new TH2F("fullMu_mass","fullMu_mass",500,0,5000,100,0,100);
};

ZDileptonAnalysis2017::ZDileptonAnalysis2017(const edm::ParameterSet& iConfig)
{
  fileName_ = iConfig.getParameter<string>("fileName");
  isMC_ = iConfig.getParameter<bool>("isMC");
  genParticleTag_ = consumes< edm::View<reco::GenParticle> >( iConfig.getParameter<edm::InputTag>("genParticleTag") );
  muTag_ = consumes< edm::View<PileupSummaryInfo> >( iConfig.getParameter<edm::InputTag>("muTag") );
  muonTag_ = consumes< edm::View<pat::Muon> >( iConfig.getParameter<edm::InputTag>("muonTag") );
  electronTag_ = consumes< edm::View<pat::Electron> >( iConfig.getParameter<edm::InputTag>("electronTag") );
  eleVetoIdMapToken_ = consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("eleVetoIdMap") );
}

ZDileptonAnalysis2017::~ZDileptonAnalysis2017()
{
 

}

// ------------ method called once each job just before starting event loop  ------------
void 
ZDileptonAnalysis2017::beginJob()
{
  root_file = new TFile(fileName_, "RECREATE");
  tree = new TTree("T", "Analysis Tree");
  totalEvts.assign(1, 0);

  tree->Branch("run", &run, "run/I");
  tree->Branch("lumi", &lumi, "lumi/I");
  tree->Branch("bx", &bx, "bx/I");
  tree->Branch("event", &event, "event/l");
  tree->Branch("rho", &rho, "rho/F");
  tree->Branch("mu", &mu, "mu/F");
  tree->Branch("nPV", &nPV, "nPV/I");
  tree->Branch("nPVall", &nPVall, "nPVall/I");
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

  double mass_ttbar = 0;
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

      if (id==6 && 20<=status && status<30) t_p4 = p.p4();
      else if (id==-6 && 20<=status && status<30) tbar_p4 = p.p4();

      else if (fabs(id)==6 && nDaught==2) {   //last t's
        const reco::GenParticle& daught0 = genParticles->at( p.daughterRef(0).key() );

        if ( fabs(daught0.pdgId())==5 || fabs(daught0.pdgId())==24 ) {
          if (id==6) t_pt2 = p.pt();
          else       tbar_pt2 = p.pt();
        }
      }
    }
    float t_pt = t_p4.Pt(), tbar_pt = tbar_p4.Pt();
    mass_ttbar = (t_p4+tbar_p4).M() > 5000 ? 4999 : (t_p4+tbar_p4).M();

    mass_totalEvts->Fill(mass_ttbar);
    ttbar_pt->Fill(t_pt, tbar_pt);
    ttbar_pt2->Fill(t_pt2, tbar_pt2);
    deltat_pt->Fill(t_pt-t_pt2);
    deltaTbar_pt->Fill(tbar_pt-tbar_pt2);


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


  //edm::Handle<edm::ValueMap<vid::CutFlowResult> > Veto_id_decisions;
  //iEvent.getByToken(eleVetoIdMapToken_, Veto_id_decisions);

  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);

  vector<pair<reco::CandidatePtr, char> > leps;

  for (int i=0, n=muons->size(); i<n; i++) {
    if (!muons->at(i).isLooseMuon()) continue;
    leps.push_back( make_pair(muons->ptrAt(i),'m') );
  }
  for (int i=0, n=electrons->size(); i<n; i++) {
    const Ptr<pat::Electron> elPtr(electrons, i);
    bool isPassEleId  = (*veto_id_decisions)[elPtr];
    cout << isPassEleId << endl;
    //cout << "for eles veto " << isPassEleId << endl;
    //if (!(*Veto_id_decisions)[elPtr].cutFlowPassed()) { cout << "non veto ele " << endl; continue;}
    //leps.push_back( make_pair(electrons->ptrAt(i),'e') );
  }

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
    if (isMC_) {
      mass_totalEvts->Write();
      ttbar_pt->Write();
      ttbar_pt2->Write();
      deltat_pt->Write();
      deltaTbar_pt->Write();
    }
    root_file->WriteObject(&totalEvts, "totalEvts");


    root_file->Write();
    delete root_file;
    root_file = 0;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZDileptonAnalysis2017);
