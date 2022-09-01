#define MAXDATA 10000
#include <iostream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>

#include "computing.h"
#include "definition.h"
#include "utility.h"

using namespace std;

double mmTomRad(double r, double inPath, double zMir){
  double path=0;
  if(inPath > 1000){
    path=inPath-zMir;
  }else{
    path=inPath+zMir;
  }
  return 1000*atan(r/path);
}

void convertToRadiant(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");
  int nedge;
  double r[MAXDATA];
  bool cutPhotonFlag[MAXDATA],outerPhoton[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("outerPhoton",&outerPhoton);

  double rRad[MAXDATA];
  auto trRad=t->Branch("rRad",&rRad,"&rRad[nedge]/D");

  cout <<"Convert radius to radiant\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++){
      double path=0;
      if(outerPhoton[j]==true){
        path = run->firstPath - run->firstMirrorPosition;
      }else{
        path = run->secondPath - run->secondMirrorPosition;
      }
      rRad[j] = atan(r[j]/path);
    }
    trRad->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}

void applyFit(TH1D *h, TF1 *f,string fname, bool out){
  if(out == false){
    h->GetXaxis()->SetRangeUser(25,50);
    f = new TF1(fname.c_str(),"gaus(0)",25,50);
  }else{
    h->GetXaxis()->SetRangeUser(130,190);
    f = new TF1(fname.c_str(),"gaus(0)",130,190);
  }
  double p1= h->GetBinCenter(h->GetMaximumBin());
  f->SetParameters(5000,p1,2);
  //h->Fit(f->GetName(),"Q","",p1-1,p1+1.5);
  h->Fit(f->GetName(),"Q","",p1-2,p1+3);
  h->Draw();
}

void computeRMS(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pmt[MAXDATA];
  double nr[MAXDATA], nt[MAXDATA], spnRadius[10], spnTime[10];
  bool coincPhoton[MAXDATA],outerPhoton[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);
  t->SetBranchAddress("spnRadius",&spnRadius);
  t->SetBranchAddress("spnTime",&spnTime);

  double rmsRadius[10], rmsTime[10], rmsPhoton[10];
  double rsdRadius[MAXDATA], rsdTime[MAXDATA];
  bool goodRMS[10];
  auto trmsRadius=t->Branch("rmsRadius",&rmsRadius,"rmsRadius[10]/D");
  auto trmsPhoton=t->Branch("rmsPhoton",&rmsPhoton,"rmsPhoton[10]/I");
  auto trmsTime=t->Branch("rmsTime",&rmsTime,"rmsTime[10]/D");
  auto tgoodRMS=t->Branch("goodRMS",&goodRMS,"goodRMS[10]/O");
  auto trsdRadius=t->Branch("rsdRadius",&rsdRadius,"rsdRadius[nedge]/D");
  auto trsdTime=t->Branch("rsdTime",&rsdTime,"rsdTime[nedge]/D");

  cout <<"Computing RMS and residui.\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < 10; j++){
      rmsRadius[j]=0;
      rmsPhoton[j]=0;
      rmsTime[j]=0;
      goodRMS[j]=false;
    }
    for(int j = 0; j < nedge; j++){
      int k=0;
      if(outerPhoton[j]==true) k = 1;
      int refPMT = pmt[j]+5*k;
      int refTOT = 4+5*k;
      rsdRadius[j]=-1000;
      rsdTime[j]=-1000;
      if(coincPhoton[j]==true){
        rsdRadius[j]=nr[j]-spnRadius[refTOT];
        rsdTime[j]=nt[j]-spnTime[refTOT];
        rmsRadius[refPMT]+=pow(spnRadius[refPMT]-nr[j],2);
        rmsTime[refPMT]+=pow(spnTime[refPMT]-nt[j],2);
        rmsPhoton[refPMT]+=1;

        rmsRadius[refTOT]+=pow(rsdRadius[j],2);
        rmsTime[refTOT]+=pow(rsdTime[j],2);
        rmsPhoton[refTOT]+=1;
      }
    }
    for(int j = 0; j < 10; j++){
      if(rmsPhoton[j]>1){
        rmsRadius[j]=sqrt(rmsRadius[j]/rmsPhoton[j]);
        rmsTime[j]=sqrt(rmsTime[j]/rmsPhoton[j]);
        goodRMS[j]=true;
      }else{
        rmsRadius[j]=-10;
        rmsTime[j]=-10;
      }
    }
    trmsRadius->Fill();
    trmsPhoton->Fill();
    trmsTime->Fill();
    tgoodRMS->Fill();
    trsdRadius->Fill();
    trsdTime->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}


void computeCutSingleParticle(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");
  int nedge, pmt[MAXDATA];
  double nr[MAXDATA], nt[MAXDATA];
  bool cutPhotonFlag[MAXDATA],outerPhoton[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("outerPhoton",&outerPhoton);
  t->SetBranchAddress("cutPhotonFlag",&cutPhotonFlag);

  double cutRadius[10], cutTime[10];
  int cutPhoton[10];
  bool goodCUT[10];
  auto tcutRadius= t->Branch("cutRadius",&cutRadius,"cutRadius[10]/D");
  auto tcutTime= t->Branch("cutTime",&cutTime,"cutTime[10]/D");
  auto tcutPhoton= t->Branch("cutPhoton",&cutPhoton,"cutPhoton[10]/I");
  auto tgoodCUT=t->Branch("goodCUT",&goodCUT,"goodCUT[10]/O");

  cout <<"Computing single particles property after cut based on rms.\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < 10; j++){
      cutRadius[j]=0;
      cutTime[j]=0;
      cutPhoton[j]=0;
    }
    for(int j = 0; j < nedge; j++){
      int k=0;
      if(outerPhoton[j]==true) k = 1;
      int refPMT = pmt[j]+5*k;
      int refTOT = 4+5*k;
      if(cutPhotonFlag[j]==true){
        cutRadius[refPMT]+=nr[j];
        cutPhoton[refPMT]+=1;
        cutTime[refPMT]+=nt[j];
        cutRadius[refTOT]+=nr[j];
        cutPhoton[refTOT]+=1;
        cutTime[refTOT]+=nt[j];
      }
    }
    for(int j = 0; j < 10; j++){
      if(cutPhoton[j]!=0){
        cutRadius[j]/=cutPhoton[j];
        cutTime[j]/=cutPhoton[j];
        goodCUT[j]=true;
        //cout <<j<<" " <<cutRadius[j]<<" " <<cutTime[j] <<" "<<cutPhoton[j] <<endl;
      }else{
        cutRadius[j]=0;
        cutTime[j]=0;
      }
    }
    tcutRadius->Fill();
    tcutTime->Fill();
    tcutPhoton->Fill();
    tgoodCUT->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}

void newSingleParticle(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);
  TFile *fIn = new TFile(fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");
  int nedge, pmt[MAXDATA];
  double nr[MAXDATA], nt[MAXDATA];
  bool coincPhoton[MAXDATA],outerPhoton[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);

  double spnRadius[10], spnTime[10];
  int spnPhoton[10];
  bool goodSPN[10];

  auto tSpnRadius=t->Branch("spnRadius",&spnRadius,"spnRadius[10]/D");
  auto tSpnPhoton=t->Branch("spnPhoton",&spnPhoton,"spnPhoton[10]/I");
  auto tSpnTime=t->Branch("spnTime",&spnTime,"spnTime[10]/D");
  auto tgoodSPN=t->Branch("goodSPN",&goodSPN,"goodSPN[10]/O");

  cout <<"Computing mean quantities for each single particle\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < 10; j++){
      spnRadius[j]=0;
      spnPhoton[j]=0;
      spnTime[j]=0;
      goodSPN[j]=false;
    }
    for(int j = 0; j < nedge; j++){
      int k=0;
      if(outerPhoton[j]==true) k = 1;
      int refPMT = pmt[j]+5*k;
      int refTOT = 4+5*k;
      if(coincPhoton[j]==true){
        spnRadius[refPMT]+=nr[j];
        spnPhoton[refPMT]+=1;
        spnTime[refPMT]+=nt[j];
        spnRadius[refTOT]+=nr[j];
        spnPhoton[refTOT]+=1;
        spnTime[refTOT]+=nt[j];
      }
    }
    for(int j = 0; j < 10; j++){
      if(spnPhoton[j]!=0){
        spnRadius[j]/=spnPhoton[j];
        spnTime[j]/=spnPhoton[j];
        goodSPN[j]=true;
        //cout <<j<<" " <<spnRadius[j]<<" " <<spnTime[j] <<" "<<spnPhoton[j] <<endl;
      }else{
        spnRadius[j]=0;
        spnTime[j]=0;
      }
    }
    tSpnRadius->Fill();
    tSpnPhoton->Fill();
    tSpnTime->Fill();
    tgoodSPN->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}


void singleParticle(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pol[MAXDATA], pmt[MAXDATA];
  double x[MAXDATA],y[MAXDATA],r[MAXDATA], nt[MAXDATA];
  bool coincPhoton[MAXDATA], outerPhoton[MAXDATA];

  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);

  double spRadius[10], spTime[10];
  int spPhoton[10];
  bool goodSP[10];

  auto tSPRadius=t->Branch("spRadius",&spRadius,"spRadius[10]/D");
  auto tSPPhoton=t->Branch("spPhoton",&spPhoton,"spPhoton[10]/I");
  auto tSPTime=t->Branch("spTime",&spTime,"spTime[10]/D");
  auto tgoodSP=t->Branch("goodSP",&goodSP,"goodSP[10]/O");

  cout <<"Selecting the photon\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < 10; j++){
      spRadius[j]=0;
      spPhoton[j]=0;
      spTime[j]=0;
      goodSP[j]=false;
    }
    for(int j = 0; j < nedge; j++){
      if(coincPhoton[j]==true){
        int k=0;
        if(outerPhoton[j]==true) k = 1;
        int refPMT = pmt[j]+5*k;
        int refTOT = 4+5*k;
        spRadius[refPMT]+=r[j];
        spPhoton[refPMT]+=1;
        spTime[refPMT]+=nt[j];
        spRadius[refTOT]+=r[j];
        spPhoton[refTOT]+=1;
        spTime[refTOT]+=nt[j];
      }
    }
    for(int j = 0; j < 10; j++){
      if(spPhoton[j]!=0){
        spRadius[j]/=spPhoton[j];
        spTime[j]/=spPhoton[j];
        goodSP[j]=true;
        //cout <<j<<" " <<spRadius[j]<<" " <<spTime[j] <<" "<<spPhoton[j] <<endl;
      }else{
        spRadius[j]=0;
        spTime[j]=0;
      }
    }
    tSPRadius->Fill();
    tSPPhoton->Fill();
    tSPTime->Fill();
    tgoodSP->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}

