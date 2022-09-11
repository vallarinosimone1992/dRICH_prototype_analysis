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

  int evt, tPol[MAXDATA], tSlot[MAXDATA], tFiber[MAXDATA], tCh[MAXDATA] ,tTime[MAXDATA], board[MAXDATA], chip[MAXDATA];
  uint tNedge;
  double tTrigTime;
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA],nt[MAXDATA];
  int pmt[MAXDATA], channel[MAXDATA];
  bool trigSig[MAXDATA];

  auto evento=tout->Branch("evt",&evt,"evt/I");
  (void)evento;
  auto ttrigtime=tout->Branch("trigtime",&tTrigTime,"trigtime/D");
  (void)ttrigtime;
  auto tnedge=tout->Branch("nedge",&tNedge,"nedge/i");
  (void)tnedge;
  auto tpol=tout->Branch("pol",&tPol[0],"pol[nedge]/I");
  (void)tpol;
  auto tslot=tout->Branch("slot",&tSlot[0],"slot[nedge]/I");
  (void)tslot;
  auto tfiber=tout->Branch("fiber",&tFiber[0],"fiber[nedge]/I");
  (void)tfiber;
  auto tch=tout->Branch("ch",&tCh[0],"ch[nedge]/I");
  (void)tch;
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
  auto tnt=tout->Branch("nt",&nt,"nt[nedge]/D");
  (void)tnt;
  auto ttrigSig=tout->Branch("trigSig",&trigSig,"trigSig[nedge]/O");

  cout <<"Taking the dRICH prototype data\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
   
    int upEdge=0;
    for(int j = 0; j < nedge; j++){
      if(pol[j]==0)upEdge++;
    }
    for(int j = 0; j < nedge; j++){
      for(int k = 0; k < nedge; k++){
        if(j == k)continue;
        if((time[j] < time[k] && pol[j]==pol[k]) || pol[j] < pol[k] ){
          int tmpTime=time[j];
          int tmpPol=pol[j];
          int tmpCh=chM[j];
          int tmpSlot=slot[j];
          int tmpFiber=fiber[j];
          time[j]=time[k];
          pol[j]=pol[k];
          chM[j]=chM[k];
          slot[j]=slot[k];
          fiber[j]=fiber[k];
          time[k]=tmpTime;
          pol[k]=tmpPol;
          chM[k]=tmpCh;
          slot[k]=tmpSlot;
          fiber[k]=tmpFiber;
        }
      }
    }
    
    evt=(int)evtDRICH;
    tTrigTime=trigtime;
    tNedge=nedge;
    bool wrongEvent=false;
    int ntrig=0;
    double trigger=0;
    for(int j = 0; j < nedge; j++){
    	if(fiber[j]==11 && chM[j]==27 && pol[j]==0){
		ntrig++;
		trigger=time[j];
	}
    }
    if(ntrig!=1)continue;
    for(uint j = 0; j < nedge; j++){
      tPol[j]=pol[j];
      tSlot[j]=slot[j];
      tFiber[j]=fiber[j];
      tCh[j]=chM[j];
      //if(fiber[j]!=11 || chM[j]!=27)tTime[j]=time[j]-trigger + 400; //Added offset 400 to distinguish peak from trigger
      //else tTime[j]=time[j]-trigger; //No offset 400 -> Trigger hit
      chip[j]=getMarocChip(chM[j]);
      board[j]=getMarocBoard(fiber[j],runHead);
      upstreamMaroc(fiber[j], runHead);
      channel[j]= getMAPMT_ch(fiber[j],chM[j],board[j],chip[j],runHead->upstreamBoard);
      if(fiber[j]!=11 || chip[j]!=0){
        tTime[j]=time[j]-trigger + 400; //Added offset 400 to distinguish peak from trigger
        trigSig[j]=false;
      }else {
        tTime[j]=time[j]-trigger; //No offset 400 -> Trigger hit
        trigSig[j]=true;
      }
      pmt[j]=FiberToPhDet(fiber[j],&(runHead->fiberRef)[0]);
      MAPMTposition(channel[j],pmt[j],&x[j],&y[j],&radius[j]);
      nt[j]=timeCalibrationMAPMT(tTime[j],channel[j],pmt[j]);
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


  int evt, tPol[MAXDATA], tSlot[MAXDATA], tFiber[MAXDATA], tCh[MAXDATA], tTime[MAXDATA], board[MAXDATA], chip[MAXDATA];
  uint tNedge;
  double tTrigTime;
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA], nt[MAXDATA];
  int pmt[MAXDATA], channel[MAXDATA];

  auto evento=tout->Branch("evt",&evt,"evt/I");
  (void)evento;
  auto ttrigtime=tout->Branch("trigtime",&tTrigTime,"trigtime/D");
  (void)ttrigtime;
  auto tnedge=tout->Branch("nedge",&tNedge,"nedge/i");
  (void)tnedge;
  auto tpol=tout->Branch("pol",&tPol[0],"pol[nedge]/I");
  (void)tpol;
  auto tslot=tout->Branch("slot",&tSlot[0],"slot[nedge]/I");
  (void)tslot;
  auto tfiber=tout->Branch("fiber",&tFiber[0],"fiber[nedge]/I");
  (void)tfiber;
  auto tch=tout->Branch("ch",&tCh[0],"ch[nedge]/I");
  (void)tch;
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
  auto tnt=tout->Branch("nt",&nt,"nt[nedge]/D");
  (void)tnt;

  cout <<"Taking the dRICH prototype data\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);

    int upEdge=0;
    for(int j = 0; j < nedge; j++){
      if(pol[j]==0)upEdge++;
    }
    for(int j = 0; j < nedge; j++){
      for(int k = 0; k < nedge; k++){
        if(j == k)continue;
        if((time[j] < time[k] && pol[j]==pol[k]) || pol[j] < pol[k] ){
          int tmpTime=time[j];
          int tmpPol=pol[j];
          int tmpCh=chM[j];
          int tmpSlot=slot[j];
          int tmpFiber=fiber[j];
          time[j]=time[k];
          pol[j]=pol[k];
          chM[j]=chM[k];
          slot[j]=slot[k];
          fiber[j]=fiber[k];
          time[k]=tmpTime;
          pol[k]=tmpPol;
          chM[k]=tmpCh;
          slot[k]=tmpSlot;
          fiber[k]=tmpFiber;
        }
      }
    }

    evt=(int)evtDRICH;
    tTrigTime=trigtime;
    tNedge=nedge;
    bool wrongEvent=false;
    int ntrig=0;
    double trigger=0;
    for(int j = 0; j < nedge; j++){
      if(fiber[j]==11 && chM[j]==27 && pol[j]==0){
        ntrig++;
        trigger=time[j];
      }
    }
    for(uint j = 0; j < nedge; j++){
      tPol[j]=pol[j];
      tSlot[j]=slot[j];
      tFiber[j]=fiber[j];
      tCh[j]=chM[j];
      if(fiber[j]!=11 || chM[j]!=27)tTime[j]=time[j]-trigger + 400; //Added offset 400 to distinguish peak from trigger
      else tTime[j]=time[j]-trigger; //No offset 400 -> Trigger hit
                                     //tTime[j]=time[j];
      board[j]=getMarocBoard(fiber[j],runHead);
      chip[j]=getMarocChip(chM[j]);
      upstreamMaroc(fiber[j], runHead);
      channel[j]= getMPPC_ch(fiber[j],chM[j],board[j],chip[j],runHead->upstreamBoard);
      pmt[j]=FiberToPhDet(fiber[j],&(runHead->fiberRef)[0]);
      MPPCposition(channel[j],pmt[j],&x[j],&y[j],&radius[j]);
      nt[j]=timeCalibrationMPPC(tTime[j],channel[j],pmt[j]);
    }
    if(wrongEvent==true)continue;
    tout->Fill();
  }
  printEnd();
  tout->Write();
  fOut->Close();
}


//void getSIMULATION(THeader *run);
