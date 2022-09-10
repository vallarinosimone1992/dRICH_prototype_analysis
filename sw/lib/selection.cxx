#define MAXDATA 10000
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>

#include "definition.h"
#include "utility.h"
#include "correction.h"

using namespace std;

////////////////////////////////////////////////////
void recoHit(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pmt[MAXDATA], pol[MAXDATA], slot[MAXDATA], fiber[MAXDATA], ch[MAXDATA];
  double x[MAXDATA], y[MAXDATA], nt[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("slot",&slot);
  t->SetBranchAddress("fiber",&fiber);
  t->SetBranchAddress("ch",&ch);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("nt",&nt);

  bool goodHit[MAXDATA];
  auto tgoodHit= t->Branch("goodHit",&goodHit,"goodHit[nedge]/O");
    
  cout <<"Applying selection based on time and radius RMS\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++)goodHit[j]=false;
    for(int j = 0; j < nedge; j++){
      //cout <<pol[j] <<endl;
      if(pol[j]!=0)continue;
      for(int k = 0; k < nedge; k++){
        if(pol[k]!=1)continue;
        //cout <<"Pol " <<pol[j] <<" " <<pol[k]<<endl;
        if(slot[j]!=slot[k])continue;
        if(fiber[j]!=fiber[k])continue;
        if(ch[j]!=ch[k])continue;
        //cout <<"Pmt " <<pmt[j] <<" " <<pmt[k]<<endl;
        if(nt[k] < nt[j])continue;
        //cout <<"Time " <<nt[j] <<" " <<nt[k] <<endl;
        if(nt[j] - nt[j] > run->MaxHitLength) continue;
        //cout <<"X " <<x[j] <<" " <<x[k] <<endl;
        //cout <<"Y " <<y[j] <<" " <<y[k] <<endl;
        if(x[j]==x[k] && y[j]==y[k]){
          goodHit[j]=true;
          goodHit[k]=true;
          nt[j] = timeWalkCorrection(nt[j],nt[k]);
          break;
        }
      }
      //cout <<"Nothing";
      //cin.get();
    }
    /*for(int j = 0; j < nedge; j++){
      if(goodHit[j]==true){
    //cout<<i <<" " <<goodHit[j] <<endl;
    //cin.get();
    }
    }*/
    tgoodHit->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}

////////////////////////////////////////////////////



void rmsCutSelection(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pmt[MAXDATA];
  double nr[MAXDATA], nt[MAXDATA];
  bool coincPhoton[MAXDATA],outerPhoton[MAXDATA];
  double rsdRadius[MAXDATA], rsdTime[MAXDATA];
  bool goodRMS[10];
  bool goodHit[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);
  t->SetBranchAddress("rsdRadius",&rsdRadius);
  t->SetBranchAddress("rsdTime",&rsdTime);
  t->SetBranchAddress("goodRMS",&goodRMS);
  t->SetBranchAddress("goodHit",&goodHit);
  t->SetBranchAddress("goodHit",&goodHit);

  bool cutPhotonFlag[MAXDATA];
  auto tcutPhotonFlag= t->Branch("cutPhotonFlag",&cutPhotonFlag,"cutPhotonFlag[nedge]/O");

  cout <<"Applying selection based on time and radius RMS\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++){
      cutPhotonFlag[j]=false;
      if(goodHit[j]==false)continue;
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
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}


////////////////////////////////////////////////////

void selectPhotons(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pol[MAXDATA];
  double x[MAXDATA],y[MAXDATA],r[MAXDATA],nt[MAXDATA];
  bool goodHit[MAXDATA];

  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("goodHit",&goodHit);

  bool coincPhoton[MAXDATA], outerPhoton[MAXDATA];
  auto tCoincPhoton= t->Branch("coincPhoton",&coincPhoton,"coincPhoton[nedge]/O");
  auto tOuterPhoton= t->Branch("outerPhoton",&outerPhoton,"outerPhoton[nedge]/O");

  cout <<"Selecting the photon in the time coincidence window and dividing rings\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++){
      coincPhoton[j]=false;
      outerPhoton[j]=false;
      if(goodHit[j]==false)continue;
      if(pol[j]==0 && nt[j] > run->timeMin && nt[j] < run->timeMax) coincPhoton[j]=true;
      //if(r[j] > run->geoCut)outerPhoton[j]=true;
      if(r[j] > run->radCut)outerPhoton[j]=true;
      //Attention, if you refine the cut you should check the conversion to milliradiant.
    }
    tCoincPhoton->Fill();
    tOuterPhoton->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}


/////////////////////////////////////////////////////

void findTimeCoincidence(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);
  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pol[MAXDATA];
  double nt[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("nt",&nt);

  TH1D *h = new TH1D("h","h",1000,0,1000);

  cout <<"Finding the time-coincidence window\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++){
      if(pol[j]==0) h->Fill(nt[j]);
    }
  }
  printEnd();
  TF1 *f0 = new TF1("f0","pol0",0,1000);
  h->Fit(f0,"Q","",800,1000);
  TF1 *f = new TF1("f","gaus(0)+pol0(3)",0,1000);
  f->SetParameters(10000,h->GetBinCenter(h->GetMaximumBin()),5,f0->GetParameter(0));
  f->SetParLimits(2,0.5,30);
  f->FixParameter(3,f0->GetParameter(0));
  double timeFitMin=h->GetBinCenter(h->GetMaximumBin())-20;
  double timeFitMax=h->GetBinCenter(h->GetMaximumBin())+20;
  h->Fit("f","Q","",timeFitMin,timeFitMax);
  //Current definition of the coincidence time window. Does it work?
  run->timeMin=f->GetParameter(1)-3*f->GetParameter(2); 
  run->timeMax=f->GetParameter(1)+3*f->GetParameter(2);
  fIn->Close();

}
