#include <iostream>
#include <iomanip>
#include <vector>
using std::cout;
using std::endl;

#include "TROOT.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

//#include "/afs/cern.ch/user/w/wulsin/root/tdrstyle.C"

#ifndef __CINT__
#endif

// ----------------------------------------------                                                                                        
// -- Global variables                         --                                                                                        
// -- Set these when running randomCommands(). --                                                                                        
// ----------------------------------------------                                                                                        

//TString inDir_ = "sigStudies/inputPlots/";
//TString outdir = "sigStudies/outputPlots/";
TCanvas* can;

// -------------------------                                                                                                             
// -- Function prototypes --                                                                                                             
// -------------------------  
void makeEff();

// -------------------------                                                                                                             
// -- Functions           --                                                                                                             
// -------------------------        

void makeEff() {



  cout << "Beginning makeEff()." << endl;

  can = new TCanvas("can", "can");
  can->Print("makeEff.ps[");

  can->SetRightMargin(0.13);
  can->SetLogx(0);
  gStyle->SetOptStat(00000000000);


  //TFile* fPi = new TFile("outputPi.root", "READ");
  TFile* fPi = new TFile("outTau.root", "READ");
  TFile* fElec = new TFile("outElec1.root", "READ");
  TFile* fMu = new TFile("outMu.root", "READ");
  
  TFile* fSel = new TFile("/afs/cern.ch/work/j/jbrinson/public/disappTrks/AnTest/CMSSW_6_1_2/src/DisappTrks/StandardAnalysis/test/condor/bkgdEstimateIs5/WjetsHighPt.root", "READ");
  
  //hists for pions
  TH1D * ptAllPi;
  TH1D * ptAllTau;

  TH1D * ptTauSel;
  TH1D * ptPiSel;

  TH1D * piPtScale;

  TH1D * pt50;
  TH1D * pt50Eff;

  TH1D * calo10Eff;
  TH1D * calo10;

  TH1D * productEff;
  TH1D * yield;

  //hists for electrons
  TH1D * ptAllElec;
  TH1D * pt50Elec;
  TH1D * pt50EffElec;
  TH1D * productElec;
  TH1D * yieldSys;

  //hists for muons
  TH1D * ptAllMu;
  TH1D * pt50Mu;
  TH1D * pt50EffMu;
  TH1D * productMu;
  TH1D * yieldSysMu;

  fPi   ->GetObject("demo/genPartPtAllPion", ptAllPi);
  fPi   ->GetObject("demo/genPartPtAllTau" , ptAllTau);
  fPi   ->GetObject("demo/genPartPtGT50"   , pt50);
  fPi   ->GetObject("demo/genPartPtLT10"   , calo10);

  fElec ->GetObject("demo/genPartPtAllPion", ptAllElec);
  fElec ->GetObject("demo/genPartPtGT50"   , pt50Elec);

  fMu ->GetObject("demo/genPartPtAllPion", ptAllMu);
  fMu ->GetObject("demo/genPartPtGT50"   , pt50Mu);

  fSel  ->GetObject("OSUAnalysis/FullSelectionTau/trackPtTrue", ptTauSel);

  //get errors
  ptAllPi   ->Sumw2();
  ptAllTau  ->Sumw2();
  pt50      ->Sumw2();
  calo10    ->Sumw2();

  ptAllElec ->Sumw2();
  pt50Elec  ->Sumw2();

  ptTauSel  ->Sumw2();

  //rebin
  ptAllPi->Rebin(20);
  ptAllTau->Rebin(20);
  pt50->Rebin(20);
  calo10->Rebin(20);

  ptAllElec->Rebin(20);
  pt50Elec->Rebin(20);

  ptAllMu->Rebin(20);
  pt50Mu->Rebin(20);

  ptTauSel->Rebin(20);

  //Calculte eff for reco pT > 50 (pions) 
  pt50Eff   = (TH1D*) pt50->Clone("pt50Eff");
  pt50Eff->SetTitle(";p^{true}_{T};P(p^{reco}_{T} > 50 GeV)");

  pt50Eff ->Sumw2();
  pt50Eff ->Divide(ptAllPi);

  //Calculte eff for reco pT > 50 (electron) 
  pt50EffElec   = (TH1D*) pt50Elec->Clone("pt50EffElec");
  pt50EffElec->SetTitle(";electron p^{true}_{T};P(electron p^{reco}_{T} > 50 GeV)");

  pt50EffElec ->Sumw2();
  pt50EffElec ->Divide(ptAllElec);

  //Calculte eff for reco pT > 50 (mu)                                                                                                                                               
  pt50EffMu   = (TH1D*) pt50Mu->Clone("pt50EffMu");
  pt50EffMu->SetTitle(";#mu p^{true}_{T};P(#mu p^{reco}_{T} > 50 GeV)");

  pt50EffMu ->Sumw2();
  pt50EffMu ->Divide(ptAllMu);


  //Calculate eff for Ecalo < 10 (pions) 
  calo10Eff = (TH1D*) calo10->Clone("calo10Eff");
  calo10Eff ->SetTitle(";p^{true}_{T};P(E_{calo}^{#DeltaR<0.5} < 10 GeV)");

  calo10Eff ->Sumw2();
  calo10Eff ->Divide(ptAllPi);

  //Calculate product eff (pions)
  productEff = (TH1D*) calo10Eff->Clone("productEff");
  productEff->SetTitle(";p^{true}_{T};P(E_{calo}^{#DeltaR<0.5} < 10 GeV) * P(p^{reco}_{T} > 50 GeV)");

  productEff->Sumw2();
  productEff->Multiply(pt50Eff);

  //Calculate product eff (electrons)
  productEffElec = (TH1D*) calo10Eff->Clone("productEffElec");
  productEffElec->SetTitle(";p^{true}_{T};P(E_{calo}^{#DeltaR<0.5} < 10 GeV) * P(electron p^{reco}_{T} > 50 GeV)");

  productEffElec->Sumw2();
  productEffElec->Multiply(pt50EffElec);

  //Calculate product eff (mu)                                                                                                                                                      
  productEffMu = (TH1D*) calo10Eff->Clone("productEffMu");
  productEffMu->SetTitle(";p^{true}_{T};P(E_{calo}^{#DeltaR<0.5} < 10 GeV) * P(#mu p^{reco}_{T} > 50 GeV)");

  productEffMu->Sumw2();
  productEffMu->Multiply(pt50EffMu);



  //reweight: pi pT/tau pT
  piPtScale = (TH1D*) ptAllPi->Clone("piPtScale");
  piPtScale->SetTitle(";#pi p_{T}^{true}/#tau p_{T}^{true};");

  piPtScale->Sumw2();
  piPtScale->Divide(ptAllTau);


  //get pi pT after our selection
  ptPiSel   = (TH1D*) ptTauSel->Clone("ptPiSel");
  ptPiSel->SetTitle(";#pi p{{T}^{true} after full selection;");

  ptPiSel->Sumw2();
  ptPiSel->Multiply(piPtScale);

  //Calculate yield for pions
  yield = (TH1D*) productEff->Clone("yield");
  yield->SetTitle(";p_{T}^{true};N_{W#rightarrow#tau#nu}");
  yield   ->Sumw2();
  //  yield   ->Multiply(ptPiSel);
  yield   ->Multiply(ptTauSel);


  //Calculate yield for electrons
  yieldSys = (TH1D*) productEffElec->Clone("yieldSys");
  yieldSys->SetTitle(";p_{T}^{true};N_{W#rightarrow#tau#nu} using P(electron p^{reco}_{T} > 50 GeV)");
  yieldSys   ->Sumw2();
  yieldSys->Multiply(ptPiSel);

  //Calculate yield for muons                                                                                                                                                          
  yieldSysMu = (TH1D*) productEffMu->Clone("yieldSysMu");
  yieldSysMu->SetTitle(";p_{T}^{true};N_{W#rightarrow#tau#nu} using P(#mu p^{reco}_{T} > 50 GeV)");
  yieldSysMu   ->Sumw2();
  yieldSysMu->Multiply(ptPiSel);


  std::cout << "The integral of the histogram for pions: " << yield->Integral() << std::endl;
  std::cout << "The integral of the histogram for electrons: " << yieldSys->Integral() << std::endl;
  std::cout << "The integral of the histogram for muons: " << yieldSysMu->Integral() << std::endl;
  double error = 0;
  double errorSys = 0;
  double errorSysMu = 0;
  std::cout << "The error of the histogram for pions: " << yield->IntegralAndError(1, yield->GetNbinsX(), error ) << std::endl;
  std::cout << "The error of the histogram is for electrons: " << yieldSys->IntegralAndError(1, yield->GetNbinsX(), errorSys ) << std::endl;
  std::cout << "The error of the histogram is for muons: " << yieldSysMu->IntegralAndError(1, yield->GetNbinsX(), errorSysMu ) << std::endl;

std::cout << "The error of the histogram for pions: " << error << std::endl;
std::cout << "The error of the histogram is for electrons: " << errorSys << std::endl;
std::cout << "The error of the histogram is for muons: " << errorSysMu << std::endl;


//yield->Multiply(ptPiSel);
//yieldSys->Multiply(ptPiSel);



  can->SetLogy(1);


  piPtScale->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  ptPiSel->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();


  pt50Eff->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  pt50EffElec->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  pt50EffMu->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();


  
  calo10Eff->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  productEff->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  productEffElec->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  productEffMu->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  yield->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  yieldSys->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  yieldSysMu->Draw("e1");
  can->Print("makeEff.ps");
  can->Clear();

  // Close the output file
  can->Print("makeEff.ps]");
  system("ps2pdf makeEff.ps");
  //  system("mv makeEff.pdf " + outdir);

  cout << "Wrote plots to: makeEff.pdf" << endl;


}

