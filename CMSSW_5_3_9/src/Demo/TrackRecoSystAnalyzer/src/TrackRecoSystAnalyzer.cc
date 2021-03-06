// -*- C++ -*-
//
// Package:    TrackRecoSystAnalyzer
// Class:      TrackRecoSystAnalyzer
// 
/**\class TrackRecoSystAnalyzer TrackRecoSystAnalyzer.cc Demo/TrackRecoSystAnalyzer/src/TrackRecoSystAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jessica Danielle Brinson,42 2-032,+41227662377,
//         Created:  Thu Nov  7 10:40:08 CET 2013
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <math.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <map>
#include <set>
#include <vector>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
//
// class declaration
//

class TrackRecoSystAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrackRecoSystAnalyzer(const edm::ParameterSet&);
      ~TrackRecoSystAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:


      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  edm::Handle<reco::GenParticleCollection> genParts;
  //  edm::Handle<reco::TrackCollection> smuons;
  edm::Handle<double> rhokt6CaloJetsHandle_;
  
  double getTrkIdMother (const reco::Track* track1, const reco::GenParticleCollection* genParts);
  int isMatchedMuonFromTau (const reco::Track* track1, const reco::GenParticleCollection* genParts);
  int isMatchedMuon (const reco::Track* track1, const reco::GenParticleCollection* genParts);
      // ----------member data ---------------------------
  edm::InputTag genPartTags_; //used to select what tracks to read from configuration file
  edm::InputTag muonTags_; //used to select what tracks to read from configuration file

  edm::Service<TFileService> fs;

  //  TFile * file;
  TH1D * hNumGlobalMu;
  TH1D * hNumSAMu;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackRecoSystAnalyzer::TrackRecoSystAnalyzer(const edm::ParameterSet& iConfig)
:
  genPartTags_(iConfig.getUntrackedParameter<edm::InputTag>("genParts")),
  muonTags_(iConfig.getUntrackedParameter<edm::InputTag>("muons"))
{
   //now do what ever initialization is needed

}


TrackRecoSystAnalyzer::~TrackRecoSystAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackRecoSystAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   using reco::TrackCollection;
   using reco::MuonCollection;
   using reco::GenParticleCollection;

   Handle<MuonCollection> muons;
      
   iEvent.getByLabel(muonTags_,muons);
   iEvent.getByLabel(genPartTags_,genParts);
   
   int numGlobalMuons = 0;
   int numSAMuons = 0;
   for(MuonCollection::const_iterator itMu = muons->begin(); itMu != muons->end();++itMu){
     if(itMu->pt() < 20) continue;
     if(itMu->isGlobalMuon() == 1) {
     std::cout << "The number of global muons = " << hNumGlobalMu->Integral() << std::endl;
     numGlobalMuons++;
     hNumGlobalMu->Fill(numGlobalMuons);
     }

     if(itMu->isStandAloneMuon() == 1) {
       std::cout << "The number of stand alone  muons = " << hNumSAMu->Integral() << std::endl;
       numSAMuons++;
       hNumSAMu->Fill(numSAMuons);
     }
   
}
   
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}





int
TrackRecoSystAnalyzer::isMatchedMuon (const reco::Track* track1, const reco::GenParticleCollection* genParts){
  //TrackRecoSystAnalyzer::isMatchedMuon (const reco::Muon* track1, const reco::GenParticleCollection* genParts){
  double value = -99;
  double genDeltaRLowest = 999;


  for(reco::GenParticleCollection::const_iterator itGenPart = genParts->begin();
      itGenPart != genParts->end();
      ++itGenPart) {
    if( fabs(itGenPart->pdgId()) == 15 ) continue;
    if( fabs(itGenPart->pdgId()) != 13 ) continue;

    double genDeltaRtemp = deltaR(itGenPart->eta(), itGenPart->phi(),track1->eta(), track1->phi());
        if (genDeltaRtemp < genDeltaRLowest) {
          genDeltaRLowest = genDeltaRtemp;
          if (genDeltaRLowest < 0.15) {   // Only consider it truth-matched if DeltaR<0.15.                                                              
            value = 1;
          }
          else {value = 0;}
        }
  }


  return value;

}

int
//TrackRecoSystAnalyzer::isMatchedMuonFromTau (const reco::Muon* track1, const reco::GenParticleCollection* genParts){
TrackRecoSystAnalyzer::isMatchedMuonFromTau (const reco::Track* track1, const reco::GenParticleCollection* genParts){
  double value = -99;
  double genDeltaRLowest = 999;
for(reco::GenParticleCollection::const_iterator itGenPart = genParts->begin();
    itGenPart != genParts->end();
    ++itGenPart) {
  if( fabs(itGenPart->pdgId()) != 15 ) continue;
  //  if (itGenPart->numberOfDaughters() > 2) continue;
  for (uint nDau = 0; nDau < itGenPart->numberOfDaughters(); nDau ++ ) {
    if (fabs(itGenPart->daughter(nDau)->pdgId() ) == 13) {
      double genDeltaRtemp = deltaR(itGenPart->eta(), itGenPart->phi(),track1->eta(), track1->phi());
      if (genDeltaRtemp < genDeltaRLowest) {
	genDeltaRLowest = genDeltaRtemp;
	if (genDeltaRLowest < 0.15) {   // Only consider it truth-matched if DeltaR<0.15.                                                                                        
	  value = 1;
	}
	else {value = 0;}
      }
    } //ends if for mu                                                                                                                                                         
  }   //ends loop over daughters                                                                                                                                                 

 }


return value;

}

double
TrackRecoSystAnalyzer::getTrkIdMother (const reco::Track* track1, const reco::GenParticleCollection* genParts){
  double value = -99;
  double genDeltaRLowest = 999;


  for(reco::GenParticleCollection::const_iterator itGenPart = genParts->begin();
      itGenPart != genParts->end();
      ++itGenPart) {
    double genDeltaRtemp = deltaR(itGenPart->eta(), itGenPart->phi(),track1->eta(), track1->phi());
    if (genDeltaRtemp < genDeltaRLowest) {
      genDeltaRLowest = genDeltaRtemp;
      if (genDeltaRLowest < 0.15) {
        if (itGenPart->numberOfMothers() > 0) {// Only consider it truth-matched if DeltaR<0.15                                                                                    
          double idMother = itGenPart->mother(0)->pdgId();
          value = idMother;
        }
      }
    }
  }

  return value;

}




// ------------ method called once each job just before starting event loop  ------------
void 
TrackRecoSystAnalyzer::beginJob()
{

  hNumGlobalMu                 = fs->make<TH1D>("hNumGlobalMu","hNumGlobalMu",3,-0.5,1.5);
  hNumSAMu                     = fs->make<TH1D>("hNumSAMu","hNumSAMu",3,-0.5,1.5);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackRecoSystAnalyzer::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
void 
TrackRecoSystAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TrackRecoSystAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TrackRecoSystAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TrackRecoSystAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackRecoSystAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

 //Specify that only 'tracks' is allowed
 //To use, remove the default given above and uncomment below
 //ParameterSetDescription desc;
 //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
 //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackRecoSystAnalyzer);
