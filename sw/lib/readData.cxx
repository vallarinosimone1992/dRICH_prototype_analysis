#define MAXDATA 10000

#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>

#include <TTree.h>
#include <TMath.h>
#include <TSystem.h>
#include <TFile.h>

#include "fillMAPS.h"
#include "photoDetPosition.h"
#include "definition.h"
#include "getChannel.h"
#include "utility.h"

void getMAPMT(THeader *runHead){
  if(runHead->sensor!="MAPMT"){
    cout <<"[ERROR] Inconsistency between the selected data type and logbook [SENSOR]\n";
    exit(EXIT_FAILURE);
  }
  //Take the sensor and time calibration data.
  getMaps();

  TString fNamedRICH=Form("%s/DATA/dRICH_DATA/run_%04d.root",runHead->suite.c_str(),runHead->runNum);
  TString fOutName=Form("%s/processed_data/firstStepData/run_%04d_processed.root",runHead->suite.c_str(),runHead->runNum);
  TFile *fdRICH = new TFile(fNamedRICH,"READ");
  TTree *t = (TTree*) fdRICH->Get("data");

  TFile *fOut = new TFile(fOutName,"RECREATE");
  TTree *tout = new TTree("dRICH","dRICH data");

  uint evtDRICH, nedge;
  int slot[MAXDATA], fiber[MAXDATA], chM[MAXDATA], pol[MAXDATA], time[MAXDATA];
  double trigtime;

  t->SetBranchAddress("evt",&evtDRICH);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("slot",&slot);
  t->SetBranchAddress("fiber",&fiber);
  t->SetBranchAddress("ch",&chM);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("time",&time);

  int evt, tPol[MAXDATA],tTime[MAXDATA], board[MAXDATA], chip[MAXDATA];
  uint tNedge;
  double tTrigTime;
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA],ntime[MAXDATA];
  int pmt[MAXDATA], channel[MAXDATA];

  auto evento=tout->Branch("evt",&evt,"evt/I");
  (void)evento;
  auto ttrigtime=tout->Branch("trigtime",&tTrigTime,"trigtime/D");
  (void)ttrigtime;
  auto tnedge=tout->Branch("nedge",&tNedge,"nedge/i");
  (void)tnedge;
  auto tpol=tout->Branch("pol",&tPol[0],"pol[nedge]/I");
  (void)tpol;
  auto tboard=tout->Branch("board",&board[0],"board[nedge]/I");
  (void)tboard;
  auto tchip=tout->Branch("chip",&chip[0],"chip[nedge]/I");
  (void)tchip;
  auto ttime=tout->Branch("time",&tTime,"time[nedge]/I");
  (void)ttime;
  auto tpmt=tout->Branch("pmt",&pmt,"pmt[nedge]/I");
  (void)tpmt;
  auto tchannel=tout->Branch("marocCh",&channel,"marocCh[nedge]/I");
  (void)tchannel;
  auto tx=tout->Branch("x",&x,"x[nedge]/D");
  (void)tx;
  auto ty=tout->Branch("y",&y,"y[nedge]/D");
  (void)ty;
  auto tradius=tout->Branch("radius",&radius,"radius[nedge]/D");
  (void)tradius;
  auto tntime=tout->Branch("ntime",&ntime,"ntime[nedge]/D");
  (void)tntime;

  cout <<"Taking the dRICH prototype data\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    evt=(int)evtDRICH;
    tTrigTime=trigtime;
    tNedge=nedge;
    bool wrongEvent=false;
    for(uint j = 0; j < nedge; j++){
      tPol[j]=pol[j];
      tTime[j]=time[j];
      board[j]=getMarocBoard(fiber[j],runHead);
      chip[j]=getMarocChip(chM[j]);
      upstreamMaroc(fiber[j], runHead);
      channel[j]= getMAPMT_ch(fiber[j],chM[j],board[j],chip[j],runHead->upstreamBoard);
      pmt[j]=FiberToPhDet(fiber[j],&(runHead->fiberRef)[0]);
      MAPMTposition(channel[j],pmt[j],&x[j],&y[j],&radius[j]);
      ntime[j]=timeCalibrationMAPMT(time[j],channel[j],pmt[j]);
    }
    if(wrongEvent==true)continue;
    tout->Fill();
  }
  printEnd();
  tout->Write();
  fOut->Close();
}


// ####################################################### //

void getMPPC(THeader *runHead){
  if(runHead->sensor!="MPPC"){
    cout <<"[ERROR] Inconsistency between the selected data type and logbook [SENSOR]\n";
    exit(EXIT_FAILURE);
  }
  //Take the sensor and time calibration data.
  getMaps();
  
  TString fNamedRICH=Form("%s/DATA/dRICH_DATA/run_%04d.root",runHead->suite.c_str(),runHead->runNum);
  TString fOutName=Form("%s/processed_data/firstStepData/run_%04d_processed.root",runHead->suite.c_str(),runHead->runNum);
  TFile *fdRICH = new TFile(fNamedRICH,"READ");
  TTree *t = (TTree*) fdRICH->Get("data");

  TFile *fOut = new TFile(fOutName,"RECREATE");
  TTree *tout = new TTree("dRICH","dRICH data");

  uint evtDRICH, nedge;
  int slot[MAXDATA], fiber[MAXDATA], chM[MAXDATA], pol[MAXDATA], time[MAXDATA];
  double trigtime;

  t->SetBranchAddress("evt",&evtDRICH);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("slot",&slot);
  t->SetBranchAddress("fiber",&fiber);
  t->SetBranchAddress("ch",&chM);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("time",&time);


  int evt, tPol[MAXDATA],tTime[MAXDATA], board[MAXDATA], chip[MAXDATA];
  uint tNedge;
  double tTrigTime;
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA], ntime[MAXDATA];
  int pmt[MAXDATA], channel[MAXDATA];

  auto evento=tout->Branch("evt",&evt,"evt/I");
  (void)evento;
  auto ttrigtime=tout->Branch("trigtime",&tTrigTime,"trigtime/D");
  (void)ttrigtime;
  auto tnedge=tout->Branch("nedge",&tNedge,"nedge/i");
  (void)tnedge;
  auto tpol=tout->Branch("pol",&tPol[0],"pol[nedge]/I");
  (void)tpol;
  auto tboard=tout->Branch("board",&board[0],"board[nedge]/I");
  (void)tboard;
  auto tchip=tout->Branch("chip",&chip[0],"chip[nedge]/I");
  (void)tchip;
  auto ttime=tout->Branch("time",&tTime,"time[nedge]/I");
  (void)ttime;
  auto tpmt=tout->Branch("pmt",&pmt,"pmt[nedge]/I");
  (void)tpmt;
  auto tchannel=tout->Branch("marocCh",&channel,"marocCh[nedge]/I");
  (void)tchannel;
  auto tx=tout->Branch("x",&x,"x[nedge]/D");
  (void)tx;
  auto ty=tout->Branch("y",&y,"y[nedge]/D");
  (void)ty;
  auto tradius=tout->Branch("radius",&radius,"radius[nedge]/D");
  (void)tradius;
  auto tntime=tout->Branch("ntime",&ntime,"ntime[nedge]/D");
  (void)tntime;

  cout <<"Taking the dRICH prototype data\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    evt=(int)evtDRICH;
    tTrigTime=trigtime;
    tNedge=nedge;
    bool wrongEvent=false;
    for(uint j = 0; j < nedge; j++){
      tPol[j]=pol[j];
      tTime[j]=time[j];
      board[j]=getMarocBoard(fiber[j],runHead);
      chip[j]=getMarocChip(chM[j]);
      upstreamMaroc(fiber[j], runHead);
      channel[j]= getMPPC_ch(fiber[j],chM[j],board[j],chip[j],runHead->upstreamBoard);
      pmt[j]=FiberToPhDet(fiber[j],&(runHead->fiberRef)[0]);
      MPPCposition(channel[j],pmt[j],&x[j],&y[j],&radius[j]);
      ntime[j]=timeCalibrationMPPC(time[j],channel[j],pmt[j]);
    }
    if(wrongEvent==true)continue;
    tout->Fill();
  }
  printEnd();
  tout->Write();
  fOut->Close();
}


//void getSIMULATION(THeader *run);
