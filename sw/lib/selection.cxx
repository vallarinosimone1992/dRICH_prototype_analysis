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


void selectPhotons(THeader *run){
  int runDRICH=run->runNum;
  //Find DRICH_SUITE environment variable
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
    cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
    exit(EXIT_FAILURE);
  }
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&env_var[0],runDRICH);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");
  
  int nedge, pol[MAXDATA], time[MAXDATA];
  double x[MAXDATA],y[MAXDATA],r[MAXDATA];

  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("time",&time);
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
      if(pol[j]==0 && time[j] > run->timeMin && time[j] < run->timeMax){
        coincPhoton[j]=true;
        if(r[j] > run->geoCut)outerPhoton[j]=true;
      }
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

  int nedge, pol[MAXDATA], time[MAXDATA];

  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("time",&time);

  cout <<"Finding the time-coincidence window\n";

  TH1D *h = new TH1D("h","h",1000,0,1000);

  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++){
      if(pol==0) h->Fill(time[j]);
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
  run->timeMin=f->GetParameter(1)-2*f->GetParameter(2);
  run->timeMax=f->GetParameter(1)+3*f->GetParameter(2);

}
