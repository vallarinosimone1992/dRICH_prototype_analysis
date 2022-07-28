#define MAXDATA 10000
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>

#include "computing.h"
#include "definition.h"

using namespace std;


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

  for(int i = 0; i < t->GetEntries(); i++){
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
  t->Write("",TObject::kOverwrite);
  fIn->Close();
  cout <<endl;
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
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}

