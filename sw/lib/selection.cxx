#define MAXDATA 10000
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>

#include "definition.h"

using namespace std;


void rmsCutSelection(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pmt[MAXDATA];
  double nr[MAXDATA], nt[MAXDATA];
  bool coincPhoton[MAXDATA],outerPhoton[MAXDATA];
  double rsdRadius[MAXDATA], rsdTime[MAXDATA];
  bool goodRMS[10];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);
  t->SetBranchAddress("rsdRadius",&rsdRadius);
  t->SetBranchAddress("rsdTime",&rsdTime);
  t->SetBranchAddress("goodRMS",&goodRMS);

  bool cutPhotonFlag[MAXDATA];
  auto tcutPhotonFlag= t->Branch("cutPhotonFlag",&cutPhotonFlag,"cutPhotonFlag[nedge]/O");
    
  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++){
      cutPhotonFlag[j]=false;
      int k=0;
      double cmpRadius = run->cutRadiusInRMS;
      double cmpTime = run->cutTimeInRMS;
      if(outerPhoton[j]==true) {
        k = 1;
        cmpRadius = run->cutRadiusOutRMS;
        cmpTime = run->cutTimeOutRMS;
      }
      int refTOT = 4+5*k;
      if(coincPhoton[j]==true && goodRMS[refTOT]==true){
        if(abs(rsdRadius[j]) < cmpRadius &&  abs(rsdTime[j]) < cmpTime)cutPhotonFlag[j] = true; 
      }
    } 
    tcutPhotonFlag->Fill();
    //cin.get();
  }  
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}


void selectPhotons(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pol[MAXDATA];
  double x[MAXDATA],y[MAXDATA],r[MAXDATA],nt[MAXDATA];

  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("r",&r);

  bool coincPhoton[MAXDATA], outerPhoton[MAXDATA];
  auto tCoincPhoton= t->Branch("coincPhoton",&coincPhoton,"coincPhoton[nedge]/O");
  auto tOuterPhoton= t->Branch("outerPhoton",&outerPhoton,"outerPhoton[nedge]/O");


  cout <<"Selecting the photon\n";
  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++){
      coincPhoton[j]=false;
      outerPhoton[j]=false;
      if(pol[j]==0 && nt[j] > run->timeMin && nt[j] < run->timeMax){
        coincPhoton[j]=true;
      }
      if(r[j] > run->geoCut)outerPhoton[j]=true;
    }
    tCoincPhoton->Fill();
    tOuterPhoton->Fill();
  }
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}




void findTimeCoincidence(THeader *run){
  int runDRICH=run->runNum;
  //Find DRICH_SUITE environment variable
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
    cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
    exit(EXIT_FAILURE);
  }
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&env_var[0],runDRICH);
  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pol[MAXDATA], nt[MAXDATA];

  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("nt",&nt);

  cout <<"Finding the time-coincidence window\n";

  TH1D *h = new TH1D("h","h",1000,0,1000);

  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++){
      if(pol==0) h->Fill(nt[j]);
    }
  }
  TF1 *f0 = new TF1("f0","pol0",0,1000);
  h->Fit(f0,"Q","",800,1000);
  TF1 *f = new TF1("f","gaus(0)+pol0(3)",0,1000);
  f->SetParameters(10000,380,5,f0->GetParameter(0));
  f->SetParLimits(2,0.5,30);
  f->FixParameter(3,f0->GetParameter(0));
  double timeFitMin=h->GetBinCenter(h->GetMaximumBin())-20;
  double timeFitMax=h->GetBinCenter(h->GetMaximumBin())+20;
  h->Fit("f","Q","",timeFitMin,timeFitMax);
  run->timeMin=f->GetParameter(1)-3*f->GetParameter(2);
  run->timeMax=f->GetParameter(1)+3*f->GetParameter(2);
  fIn->Close();

}
