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

void getMAPMT(THeader *runHead){
  int runDRICH=runHead->runNum;
  string phDet=runHead->sensor;
  if(phDet!="MAPMT"){
    cout <<phDet.c_str() <<endl;
    cout <<"[ERROR] Inconsistency between the selected data type and logbook [SENSOR]\n";
    exit(EXIT_FAILURE);
  }
  //Take the MAPMT maps
  map<string,int> map_MAPMT1;
  map<string,int> map_MAPMT2;
  map<string,int>::iterator it_map_MAPMT1;
  map<string,int>::iterator it_map_MAPMT2;
  getMapMAPMT(&map_MAPMT1,&map_MAPMT2);

  getMaps();

  //Find DRICH_SUITE environment variable
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
    cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
    exit(EXIT_FAILURE);
  }

  TString fNamedRICH=Form("%s/DATA/dRICH_DATA/run_%04d.root",&env_var[0],runDRICH);
  TString fOutName=Form("%s/processed_data/firstStepData/run_%04d_processed.root",&env_var[0],runDRICH);
  //TString fOutName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_processed.root",runHead->suite,runDRICH);
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

  cout <<"dRICH variables setted\n";

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

  for(int i = 0; i < t->GetEntries(); i++){
    if((i)%(t->GetEntries()/10)==0)cout <<Form("\rTaking the dRICH data: %lld%% completed    ",(1+(i/(t->GetEntries()/10)))*10) <<flush;
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
  cout <<endl;
  tout->Write();
  fOut->Close();
  cout <<Form("MAPMTs data wrote in %s\n",&fOutName[0]);
}


// ####################################################### //

void getMPPC(THeader *runHead){
  cout <<"Inside\n";
  int runDRICH=runHead->runNum;
  string phDet=runHead->sensor;
  if(phDet!="MPPC"){
    cout <<"[ERROR] Inconsistency between the selected data type and logbook [SENSOR]\n";
    exit(EXIT_FAILURE);
  }
  cout <<"Sensor checked\n";
  //Take the MAPMT maps
  map<string,int> map_MPPC1;
  map<string,int> map_MPPC2;
  map<string,int>::iterator it_map_MPPC1;
  map<string,int>::iterator it_map_MPPC2;
  getMapMPPC(&map_MPPC1,&map_MPPC2);

  //Find DRICH_SUITE environment variable
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
    cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
    exit(EXIT_FAILURE);
  }
  cout <<"DRICH_SUITE trovata\n";

  TString fNamedRICH=Form("%s/DATA/dRICH_DATA/run_%04d.root",&env_var[0],runDRICH);
  TString fOutName=Form("%s/processed_data/firstStepData/run_%04d_processed.root",&env_var[0],runDRICH);
  //TString fOutName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_processed.root",runHead->suite,runDRICH);
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

  cout <<"dRICH variables setted\n";

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

  cout <<"TTree presi\n";

  for(int i = 0; i < t->GetEntries(); i++){
    if((i)%(t->GetEntries()/10)==0)cout <<Form("\rTaking the dRICH data: %lld%% completed    ",(1+(i/(t->GetEntries()/10)))*10) <<flush;
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
  cout <<endl;
  tout->Write();
  fOut->Close();
  cout <<Form("MPPCs data wrote in %s\n",&fOutName[0]);
}


//void getSIMULATION(THeader *run);