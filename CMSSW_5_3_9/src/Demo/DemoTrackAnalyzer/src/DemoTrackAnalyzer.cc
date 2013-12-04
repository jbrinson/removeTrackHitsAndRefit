// -*- C++ -*-
//
// Package:    DemoTrackAnalyzer
// Class:      DemoTrackAnalyzer
// 
/**\class DemoTrackAnalyzer DemoTrackAnalyzer.cc Demo/DemoTrackAnalyzer/src/DemoTrackAnalyzer.cc

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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
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

class DemoTrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DemoTrackAnalyzer(const edm::ParameterSet&);
      ~DemoTrackAnalyzer();

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
  edm::Handle<reco::TrackCollection> muons;
  edm::Handle<reco::TrackCollection> smuons;
  edm::Handle<double> rhokt6CaloJetsHandle_;
  
  double getTrkPtRes (const reco::Track* track1);
  double getTrkPtPull (const reco::Track* track1);
  double getTrkPtTrue (const reco::Track* track1, const reco::GenParticleCollection* genParts);
  double getTrkPtTrueMother (const reco::Track* track1, const reco::GenParticleCollection* genParts);
  double getTrkIdMother (const reco::Track* track1, const reco::GenParticleCollection* genParts);
  int isMatchedPionFromTau (const reco::Track* track1, const reco::GenParticleCollection* genParts);
  int isMatchedTau (const reco::Track* track1, const reco::GenParticleCollection* genParts);
  int isMatchedPion (const reco::Track* track1, const reco::GenParticleCollection* genParts);
  int getTrkIsIso (const reco::Track* track1, const reco::TrackCollection* trackColl);
  double getTrkCaloTotRhoCorr(const double em, const double had);
      // ----------member data ---------------------------
  edm::InputTag trackTags_; //used to select what tracks to read from configuration file
  edm::InputTag genPartTags_; //used to select what tracks to read from configuration file
  edm::InputTag muonTags_; //used to select what tracks to read from configuration file
  edm::InputTag smuonTags_; //used to select what tracks to read from configuration file
  edm::InputTag simvtxTag_; //used to select what tracks to read from configuration file
  bool isPiPartGun_; //used to select what tracks to read from configuration file
  edm::Service<TFileService> fs;

  //  TFile * file;
  TH1D * trackPt;
  TH1D * tauMom;
  TH1D * trackPtAllPion;
  TH1D * trackPtCaloDiffAllPion;
  TH1D * trackEtaAllPion;
  TH1D * caloTotAllPion;
  TH1D * trackPtGen15;
  TH1D * trackPtGen30;
  TH1D * trackPtGen50;
  TH1D * trackPtLT;
  TH1D * trackPtCaloDiffLT;
  TH1D * trackEtaLT;
  TH1D * trackEtaGT;
  TH1D * trackPtCaloDiffGT;
  TH1D * trackPtTrue;
  TH1D * trackPtTrueLT;
  TH1D * genPt;
  TH1D * genPartPtAllPion;
  TH1D * genPartPtAllTau;
  TH1D * genPartPtGT50;
  TH1D * genPartPtGT50Eff;
  TH1D * genPartPtLT10;
  TH1D * genPartPtLT3;
  TH1D * genPartPtLT10Eff;
  TH1D * ProductEff;
  TH1D * trackPtError;
  TH1D * trackFracPtError;
  TH1D * trackPtRes;
  TH1D * trackPtPull;
  TH1D * trackNumValidHits;
  TH1D * trackNMissOut;
  TH1D * hnumPion;
  TH1D * trackCaloTot;
  TH1D * trackCaloTotPtGTFifty;
  TH2D *  hPosSimVtx;
  TH2D *  hEffVsPt;
  TH2D *  hEffVsPtTrue;
  
  TH2D *  caloTotVsPtTrue;
  TH2D *  ptVsPtTrue;
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
DemoTrackAnalyzer::DemoTrackAnalyzer(const edm::ParameterSet& iConfig)
:
  trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
  genPartTags_(iConfig.getUntrackedParameter<edm::InputTag>("genParts")),
  muonTags_(iConfig.getUntrackedParameter<edm::InputTag>("muons")),
  smuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("smuons")),
  simvtxTag_(iConfig.getUntrackedParameter< edm::InputTag >("simvtxTag")),
  isPiPartGun_(iConfig.getUntrackedParameter<bool >("isPiPartGun"))

{
   //now do what ever initialization is needed

}


DemoTrackAnalyzer::~DemoTrackAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   using reco::TrackCollection;
   using reco::GenParticleCollection;

   Handle<TrackCollection> tracks;
   Handle<SimVertexContainer> simVertexCollection;

   iEvent.getByLabel(simvtxTag_, simVertexCollection);
   iEvent.getByLabel(trackTags_,tracks);
   iEvent.getByLabel(muonTags_,muons);
   iEvent.getByLabel(smuonTags_,smuons);
   iEvent.getByLabel(genPartTags_,genParts);
   iEvent.getByLabel ("kt6CaloJets","rho", rhokt6CaloJetsHandle_);

   const SimVertexContainer simVC = *(simVertexCollection.product());

   int numPion = 0;

   int numMuons = 0;
   for(TrackCollection::const_iterator itMu = muons->begin(); itMu != muons->end();++itMu){
     numMuons++;
   }

   int numsMuons = 0;
   for(TrackCollection::const_iterator itMu = smuons->begin(); itMu != smuons->end();++itMu){
     numsMuons++;
   }


   for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end();++itTrack) {

     if(numMuons > 0) continue;
     if(numsMuons > 0) continue;

     //get associated calo energy for each trk
     std::vector<edm::Handle<CaloTowerCollection> > prods;
     try {
       iEvent.getManyByType(prods);
     } catch (...) {
       std::cout << "No CaloTowers." << std::endl;
     }
     //Access info about Calo Towers                                                                                  
     double caloEMDeltaRp5  = 0;
     double caloHadDeltaRp5 = 0;
     
     std::vector<edm::Handle<CaloTowerCollection> >::iterator i = prods.begin();
     const CaloTowerCollection& c=*(*i);
     for (CaloTowerCollection::const_iterator j=c.begin(); j!=c.end(); j++) {
       double deltaEta = fabs(itTrack->eta() - j->eta());
       double deltaPhi = fabs(fabs(fabs(itTrack->phi() - j->phi()) - ROOT::Math::Pi()) - ROOT::Math::Pi());
       double deltaR   = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
       double Eem  = j->emEnergy();
       double Ehad = j->hadEnergy();

       if (j->emEt()  < 0.2) Eem  = 0;
       if (j->hadEt() < 0.5) Ehad = 0;

       if (deltaR<0.5) {
	 caloEMDeltaRp5  += Eem;
	 caloHadDeltaRp5 += Ehad;
       }
     }


     if(isPiPartGun_) {
       if(isMatchedPion(&(*itTrack), genParts.product()) <= 0 ) continue;
     }

     if(!isPiPartGun_) {
       //       if(isMatchedPionFromTau(&(*itTrack), genParts.product()) <= 0 ) continue; 
       if(isMatchedTau(&(*itTrack), genParts.product()) <= 0 ) continue;                                                                               
     }


    if(itTrack->pt() < 5) continue;
    //    if( fabs(itTrack->eta()) < 2.1) continue;
    if(itTrack->numberOfValidHits() < 5) continue;
	//if(itTrack->numberOfValidHits() < 6) continue;
     //only consider tracks that are matched to a pi from a tau in the case of WtoTaNu                                                                   
     //     tauMom->Fill( getTrkIdMother(&(*itTrack), genParts.product()) );
     if ( fabs(getTrkIdMother(&(*itTrack), genParts.product()))  == 15) {
     genPartPtAllTau->Fill(getTrkPtTrueMother(&(*itTrack), genParts.product()));
}

     //     if (getTrkIsIso(&(*itTrack), tracks.product()) == 0) continue;
     numPion++;

     //Fill histos for all pions
     hnumPion->Fill(numPion);
     genPartPtAllPion->Fill(getTrkPtTrue(&(*itTrack), genParts.product()));
     //     trackPtCaloDiffAllPion->Fill( (itTrack->pt()) - ( caloEMDeltaRp5+caloHadDeltaRp5 ));
     trackPtCaloDiffAllPion->Fill( (itTrack->pt()) - getTrkCaloTotRhoCorr( caloEMDeltaRp5, caloHadDeltaRp5 ));
     trackPtAllPion->Fill(itTrack->pt());
     trackEtaAllPion->Fill(itTrack->eta());
     caloTotAllPion->Fill( getTrkCaloTotRhoCorr(caloEMDeltaRp5, caloHadDeltaRp5) );
     caloTotVsPtTrue->Fill( getTrkCaloTotRhoCorr(caloEMDeltaRp5, caloHadDeltaRp5), getTrkPtTrue(&(*itTrack), genParts.product()));
     ptVsPtTrue->Fill( itTrack->pt(), getTrkPtTrue(&(*itTrack), genParts.product()));

     //fill histos for pions with reco pt > 50
     if(itTrack->pt() > 50){
       genPartPtGT50->Fill(getTrkPtTrue(&(*itTrack), genParts.product()));
       trackEtaGT->Fill(itTrack->eta());
       trackPtCaloDiffGT->Fill( (itTrack->pt()) - ( getTrkCaloTotRhoCorr(caloEMDeltaRp5,caloHadDeltaRp5)));

       
     }

     //fill histos for pions with Ecalo < 10	
     if ( getTrkCaloTotRhoCorr(caloEMDeltaRp5, caloHadDeltaRp5) < 10){

       std::cout<<iEvent.id().event()<<std::endl;

       genPartPtLT10->Fill(getTrkPtTrue(&(*itTrack), genParts.product()));
       trackEtaLT->Fill(itTrack->eta());
       trackPtCaloDiffLT->Fill( (itTrack->pt()) - ( getTrkCaloTotRhoCorr(caloEMDeltaRp5,caloHadDeltaRp5) ));
       //       std::cout << "The generated pT = " << getTrkPtTrue(&(*itTrack), genParts.product()) << std::endl;
       trackPtTrueLT->Fill(getTrkPtTrue(&(*itTrack), genParts.product()));
       //std::cout << "The reconstructed pT = " << itTrack->pt() << " eta = " << itTrack->eta()<< " in event: " << iEvent.id().event() <<   std::endl;
       trackPtLT->Fill(itTrack->pt());
     }
     
     //fill reco track pt in coarse bins of gen pt
     if (getTrkPtTrue(&(*itTrack), genParts.product()) < 15 ){
       //       std::cout<<"Gen pt (<15) = " << getTrkPtTrue(&(*itTrack), genParts.product())<<std::endl; 
       trackPtGen15->Fill(itTrack->pt());
     }
     
     if (getTrkPtTrue(&(*itTrack), genParts.product()) > 15 && getTrkPtTrue(&(*itTrack), genParts.product()) < 30 ){
       //       std::cout<<"Gen pt (15to30) = " << getTrkPtTrue(&(*itTrack), genParts.product())<<std::endl; 
       trackPtGen30->Fill(itTrack->pt());
     }
     if (getTrkPtTrue(&(*itTrack), genParts.product()) > 30 && getTrkPtTrue(&(*itTrack), genParts.product()) < 50 ){
       //       std::cout<<"Gen pt (30to50) = " << getTrkPtTrue(&(*itTrack), genParts.product())<<std::endl; 
       trackPtGen50->Fill(itTrack->pt());
     }
     
     if (getTrkPtRes(&(*itTrack)) > 4){
       std::cout<<"PtRes is Greater > 4" << std::endl;
     }
     if (getTrkCaloTotRhoCorr(caloEMDeltaRp5, caloHadDeltaRp5) < 10){
       //  std::cout<< " I have a suspiciously low Ecalo. Look at event = " << iEvent.id().event() << std::endl;
     }
          
     trackCaloTot->Fill( getTrkCaloTotRhoCorr(caloEMDeltaRp5, caloHadDeltaRp5));
      trackPt->Fill(itTrack->pt());
      trackPtTrue->Fill(getTrkPtTrue(&(*itTrack), genParts.product()) );
      trackPtError->Fill(itTrack->ptError());
      trackFracPtError->Fill(itTrack->ptError()/itTrack->pt());
      trackNumValidHits->Fill(itTrack->numberOfValidHits());
      trackNMissOut->Fill(itTrack->trackerExpectedHitsOuter().numberOfHits());
      trackPtRes->Fill(getTrkPtRes(&(*itTrack)));

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





double
//DemoTrackAnalyzer::getTrkPtRes (const reco::Track* track1){
DemoTrackAnalyzer::getTrkPtRes (const reco::Track* track1){

  double ptTrue = getTrkPtTrue(track1, genParts.product());
  double PtRes = (track1->pt() - ptTrue) / ptTrue;

  //  std::cout << "pT True = " << ptTrue << std::endl;
  return PtRes;

}

double
//DemoTrackAnalyzer::getTrkPtRes (const reco::Track* track1){
DemoTrackAnalyzer::getTrkPtPull (const reco::Track* track1){

  double ptTrue = getTrkPtTrue(track1, genParts.product());
  double PtPull = (track1->pt() - ptTrue) / track1->ptError();

  return PtPull;

}


double
DemoTrackAnalyzer::getTrkPtTrue (const reco::Track* track1, const reco::GenParticleCollection* genParts){
  double value = -99;
  double genDeltaRLowest = 999;


  for(reco::GenParticleCollection::const_iterator itGenPart = genParts->begin();
      itGenPart != genParts->end();                      
      ++itGenPart) {
    genPt->Fill(itGenPart->pt());
       double genDeltaRtemp = deltaR(itGenPart->eta(), itGenPart->phi(),track1->eta(), track1->phi());
      if (genDeltaRtemp < genDeltaRLowest) {
	genDeltaRLowest = genDeltaRtemp;
	if (genDeltaRLowest < 0.15) {   // Only consider it truth-matched if DeltaR<0.15.
	  double ptTrue = itGenPart->pt();
	  value = ptTrue;
	}
    }
  }
  
  return value;
  
}

double
DemoTrackAnalyzer::getTrkPtTrueMother (const reco::Track* track1, const reco::GenParticleCollection* genParts){
  double value = -99;
  double genDeltaRLowest = 999;


  for(reco::GenParticleCollection::const_iterator itGenPart = genParts->begin();
      itGenPart != genParts->end();
      ++itGenPart) {
    genPt->Fill(itGenPart->pt());
    double genDeltaRtemp = deltaR(itGenPart->eta(), itGenPart->phi(),track1->eta(), track1->phi());
    if (genDeltaRtemp < genDeltaRLowest) {
      genDeltaRLowest = genDeltaRtemp;
      if (genDeltaRLowest < 0.15) {   
	if (itGenPart->numberOfMothers() > 0) {// Only consider it truth-matched if DeltaR<0.15
	  double ptTrueMother = itGenPart->mother(0)->pt();
	value = ptTrueMother;
	}
	
	if (itGenPart->numberOfMothers() > 1) {// Only consider it truth-matched if DeltaR<0.15                                                                
	  double ptTrueMother2 = itGenPart->mother(1)->pt();
	  value = ptTrueMother2;
        }
      }
    }
  }

  return value;

}

double
DemoTrackAnalyzer::getTrkIdMother (const reco::Track* track1, const reco::GenParticleCollection* genParts){
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





int
DemoTrackAnalyzer::getTrkIsIso (const reco::Track* track1, const reco::TrackCollection* trackColl){

  for(reco::TrackCollection::const_iterator track2 = trackColl->begin(); track2 !=trackColl->end(); track2++){
    if(track1->eta() == track2->eta() && track1->phi() == track2->phi()) continue; // Do not compare the track to itself.
    double deltaRtrk = deltaR(track1->eta(), track1->phi(), track2->eta(), track2->phi());
    if(deltaRtrk < 0.2) return 0;

  }
  return 1;

}


double
DemoTrackAnalyzer::getTrkCaloTotRhoCorr(const double em, const double had) {
  // Return the pile-up (rho) corrected isolation energy, i.e., the total calorimeter energy around the candidate track.
  double radDeltaRCone = 0.5;
  double rhoCorr_kt6CaloJets = *rhokt6CaloJetsHandle_ * ROOT::Math::Pi() * pow(radDeltaRCone, 2);  // Define effective area as pi*r^2, where r is radius of DeltaR cone.
  double rawCaloTot = em + had;
  double caloTotRhoCorrCalo = TMath::Max(0., rawCaloTot - rhoCorr_kt6CaloJets);
  return caloTotRhoCorrCalo;

}

int
DemoTrackAnalyzer::isMatchedPionFromTau (const reco::Track* track1, const reco::GenParticleCollection* genParts){
  double value = -99;
  double genDeltaRLowest = 999;


  for(reco::GenParticleCollection::const_iterator itGenPart = genParts->begin();
      itGenPart != genParts->end();                      
      ++itGenPart) {
    if( fabs(itGenPart->pdgId()) != 15 ) continue; 
    if (itGenPart->numberOfDaughters() > 2) continue;
    for (uint nDau = 0; nDau < itGenPart->numberOfDaughters(); nDau ++ ) {
      if (fabs(itGenPart->daughter(nDau)->pdgId() ) == 211) {
	double genDeltaRtemp = deltaR(itGenPart->eta(), itGenPart->phi(),track1->eta(), track1->phi());
	if (genDeltaRtemp < genDeltaRLowest) {
	  genDeltaRLowest = genDeltaRtemp;
	  if (genDeltaRLowest < 0.15) {   // Only consider it truth-matched if DeltaR<0.15.
	    value = 1;
	  }
	  else {value = 0;} 
	}
      } //ends if for pion
    }   //ends loop over daughters
    
  }
  
  
  return value;

}

int
DemoTrackAnalyzer::isMatchedTau (const reco::Track* track1, const reco::GenParticleCollection* genParts){
  double value = -99;
  double genDeltaRLowest = 999;


  for(reco::GenParticleCollection::const_iterator itGenPart = genParts->begin();
      itGenPart != genParts->end();
      ++itGenPart) {
    if( fabs(itGenPart->pdgId()) != 15 ) continue;
    if (itGenPart->numberOfDaughters() > 2) continue;
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
DemoTrackAnalyzer::isMatchedPion (const reco::Track* track1, const reco::GenParticleCollection* genParts){
  double value = -99;
  double genDeltaRLowest = 999;


  for(reco::GenParticleCollection::const_iterator itGenPart = genParts->begin();
      itGenPart != genParts->end();
      ++itGenPart) {
    if( fabs(itGenPart->pdgId()) != 211 ) continue;
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


// ------------ method called once each job just before starting event loop  ------------
void 
DemoTrackAnalyzer::beginJob()
{
  //file = new TFile("outfile.root","recreate");
  //const bool oldAddDir = TH1::AddDirectoryStatus();
  //TH1::AddDirectory(true);

trackPt                    = fs->make<TH1D>("trackPt","trackPt",100,0,150);
tauMom                     = fs->make<TH1D>("tauMom","tauMom",300,0,300);
trackPtAllPion             = fs->make<TH1D>("trackPtAllPion","trackPtAllPion",100,0,150);
trackPtCaloDiffAllPion     = fs->make<TH1D>("trackPtCaloDiffAllPion","trackPtCaloDiffAllPion",100,-75,75);
trackEtaAllPion            = fs->make<TH1D>("trackEtaAllPion","trackEtaAllPion",100,-3,3);
caloTotAllPion             = fs->make<TH1D>("caloTotAllPion","caloTotAllPion",100,0,75);
trackPtGen15               = fs->make<TH1D>("trackPtGen15","trackPtGen15",100,0,150);
trackPtGen30               = fs->make<TH1D>("trackPtGen30","trackPtGen30",100,0,150);
trackPtGen50               = fs->make<TH1D>("trackPtGen50","trackPtGen50",100,0,150);
trackPtLT                  = fs->make<TH1D>("trackPtLT","trackPtLT",100,0,150);
trackPtCaloDiffLT          = fs->make<TH1D>("trackPtCaloDiffLT","trackPtCaloDiffLT",100,-75,75);
trackPtCaloDiffGT          = fs->make<TH1D>("trackPtCaloDiffGT","trackPtCaloDiffGT",100,-75,75);
trackEtaLT                 = fs->make<TH1D>("trackEtaLT","trackEtaLT",100,-3,3);
trackEtaGT                 = fs->make<TH1D>("trackEtaGT","trackEtaGT",100,-3,3);
trackPtTrue                = fs->make<TH1D>("trackPtTrue","trackPtTrue",100,0,150);
trackPtTrueLT              = fs->make<TH1D>("trackPtTrueLT","trackPtTrue",1000,0,150);
trackCaloTot               = fs->make<TH1D>("trackCaloTot","trackCaloTot",100,0,50);
trackCaloTotPtGTFifty      = fs->make<TH1D>("trackCaloTotPtGTFifty","trackCaloTotPtGTFifty",100,0,50);
genPt                      = fs->make<TH1D>("genPt","genPt",100,0,100);
genPartPtAllPion           = fs->make<TH1D>("genPartPtAllPion","genPartPtAllPion",100,0,150);
genPartPtAllTau            = fs->make<TH1D>("genPartPtAllTau","genPartPtAllTau",100,0,150);
genPartPtGT50              = fs->make<TH1D>("genPartPtGT50","genPartPtGT50",100,0,150);
//genPartPtGT50Eff           = fs->make<TH1D>("genPartPtGT50Eff","genPartPtGT50Eff",100 ,0,100);
genPartPtLT10              = fs->make<TH1D>("genPartPtLT10","genPartPtLT10",100,0,150);
genPartPtLT3              = fs->make<TH1D>("genPartPtLT3","genPartPtLT3",100,0,150);
//genPartPtLT10Eff           = fs->make<TH1D>("genPartPtLT10Eff","genPartPtLT10Eff",100,0,100);
//ProductEff               = fs->make<TH1D>("genPartPtLT10Eff","genPartPtLT10Eff",50,0,100);
trackPtError               = fs->make<TH1D>("ptError","ptError",100,0,100);
trackFracPtError           = fs->make<TH1D>("fracPtError","fracPtError",100,0,1);
trackPtRes                 = fs->make<TH1D>("trackPtRes","trackPtRes",100,-3,3);
trackPtPull                = fs->make<TH1D>("trackPtPull","trackPtPull",100,-3,3);
trackNumValidHits          = fs->make<TH1D>("trackNumValidHits","trackNumValidHits",100,-0.5,14.5);
trackNMissOut              = fs->make<TH1D>("trackNMissOut","trackNMissOut",100,-0.5,14.5);
hnumPion                   = fs->make<TH1D>("hnumPion","hnumPion",3,-0.5,1.5);

hPosSimVtx                 = fs->make<TH2D>("hPosSimVtx",";z position (sim vertices);radius",100, -300, 300, 1000, 0, 120);
hEffVsPt                   = fs->make<TH2D>("hEffVsPt","reco pt;eff",100, 0, 300, 100, 0, 1);
hEffVsPtTrue               = fs->make<TH2D>("hEffVsPtTrue","true pt;eff",100, 0, 300, 100, 0, 1);

caloTotVsPtTrue            = fs->make<TH2D>("caloTotVsPtTrue",";caloTot;ptTrue ",100, 0, 75, 100, 0, 150);
ptVsPtTrue                 = fs->make<TH2D>("ptVsPtTrue",";pt;pt true",100, 0, 150, 100, 0, 150);

//TH1::AddDirectory(oldAddDir);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoTrackAnalyzer::endJob() 
{

  //TH1D* genPartPtGT50Eff = (TH1D*) genPartPtGT50->Clone("genPartPtGT50Eff");                                                                                                 
  genPartPtGT50Eff = (TH1D*) genPartPtGT50->Clone("genPartPtGT50Eff");
  //      std::cout<< "Bin content of eficiency plot before division = " << genPartPtGT50Eff->GetBinContent(27) << std::endl;                                                  
        std::cout<< "Integral of eficiency plot before division = " << genPartPtGT50Eff->Integral() << std::endl;                                                            
  genPartPtGT50Eff->Divide(genPartPtAllPion);
        std::cout<< "Integral of eficiency plot  = " << genPartPtGT50Eff->Integral() << std::endl;                                                                           
  //      std::cout<< "Bin content of eficiency plot = " << genPartPtGT50Eff->GetBinContent(27) << std::endl;                                                                  
  genPartPtLT10Eff = (TH1D*) genPartPtLT10->Clone("genPartPtLT10Eff");
  genPartPtLT10Eff->Divide(genPartPtAllPion);

  genPartPtGT50Eff->Write();
  genPartPtLT10Eff->Write();
  //file->Write();
  //file->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
DemoTrackAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DemoTrackAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DemoTrackAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DemoTrackAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoTrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(DemoTrackAnalyzer);
