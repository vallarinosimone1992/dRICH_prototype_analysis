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

#include "integrate.h"
#include "utility.h"
#include "tracking.h"
#include "getChannel.h"
#include "photoDetPosition.h"
#include "fillMAPS.h"
#include "definition.h"

using namespace std;


void noGEM_Integration(THeader *run){
  TString fNamedRICH=Form("%s/processed_data/firstStepData/run_%04d_processed.root",run->suite.c_str(),run->runNum);
  TString fOutName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  if(gSystem->AccessPathName(fNamedRICH)){
    cout <<"[ERROR] dRICH run not found\n";
    exit(EXIT_FAILURE);
  }
  //Getting the TTrees
  TFile *fdRICH = new TFile(fNamedRICH,"READ");
  if(fdRICH->IsZombie()){
    cout <<"[ERROR] fdRICH is zombie! Something was wrong\n";
    exit(EXIT_FAILURE);
  }
  TTree *t = (TTree*) fdRICH->Get("dRICH");
  t->SetTitle("dRICH and GEM data");

  uint nedge;
  int evt, board[MAXDATA], chip[MAXDATA], pol[MAXDATA], time[MAXDATA], marocCh[MAXDATA], pmt[MAXDATA];
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA], ntime[MAXDATA];
  double trigtime;

  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("board",&board);
  t->SetBranchAddress("chip",&chip);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("marocCh",&marocCh);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("radius",&radius);
  t->SetBranchAddress("ntime",&ntime);

  float gx0=0, gy0=0, gx1=0, gy1=0, gxa=0, gya=0, gxtheta=0, gytheta=0;
  
  TFile *fOut = new TFile(fOutName,"RECREATE"); 
  TTree *tout = new TTree("dRICH","dRICH integrated data");
  
  int tevt, tnedge, tpol[MAXDATA], tboard[MAXDATA], tchip[MAXDATA],ttime[MAXDATA],tpmt[MAXDATA], tchannel[MAXDATA];
  double ttrigtime, tx[MAXDATA],ty[MAXDATA],tr[MAXDATA], tnt[MAXDATA];
  bool dataGEM;

  auto tEvt=tout->Branch("evt",&tevt,"evt/I");
  auto tTrigtime=tout->Branch("trigtime",&ttrigtime,"trigtime/I");
  auto tNedge=tout->Branch("nedge",&tnedge,"nedge/I");
  auto tPol=tout->Branch("pol",&tpol,"pol[nedge]/I");
  auto tTime=tout->Branch("time",&ttime,"time[nedge]/I");
  auto tPmt=tout->Branch("pmt",&tpmt,"pmt[nedge]/I");
  auto tBoard=tout->Branch("board",&tboard,"board[nedge]/I");
  auto tChip=tout->Branch("chip",&tchip,"chip[nedge]/I");
  auto tX=tout->Branch("x",&tx,"x[nedge]/D");
  auto tY=tout->Branch("y",&ty,"y[nedge]/D");
  auto tR=tout->Branch("r",&tr,"r[nedge]/D");
  auto tNT=tout->Branch("nt",&tnt,"nt[nedge]/D");

  auto tGX0=tout->Branch("gx0",&gx0,"gx0/F");
  auto tGY0=tout->Branch("gy0",&gy0,"gy0/F");
  auto tGX1=tout->Branch("gx1",&gx1,"gx1/F");
  auto tGY1=tout->Branch("gy1",&gy1,"gy1/F");
  auto tGXA=tout->Branch("gxa",&gxa,"gxa/F");
  auto tGYA=tout->Branch("gya",&gya,"gya/F");
  auto tGXtheta=tout->Branch("gxtheta",&gxtheta,"gxtheta/F");
  auto tGYtheta=tout->Branch("gytheta",&gytheta,"gytheta/F");


  cout <<"No tracking info. I'm filling the tracking data with null values\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    dataGEM=false;
    t->GetEntry(i);
    //Copy the dRICH info
    tevt=(int)evt;
    ttrigtime=trigtime;
    tnedge=(int)nedge;
    for(int j = 0; j < nedge; j++){
      tpol[j]=pol[j];
      tboard[j]=board[j];
      tchip[j]=chip[j];
      ttime[j]=time[j];
      tpmt[j]=pmt[j];
      tx[j]=x[j];
      ty[j]=y[j];
      tr[j]=radius[j];
      tnt[j]=ntime[j];
    }
    tout->Fill();
  }
  printEnd(); 
  tout->Write();
  fOut->Close();
  fdRICH->Close();
}

void TTreeIntegration(THeader *run){
  TString fNamedRICH=Form("%s/processed_data/firstStepData/run_%04d_processed.root",run->suite.c_str(),run->runNum);
  TString fNameGEM=Form("%s/DATA/GEM_DATA/run_%04d_gem.root",run->suite.c_str(),run->runNumGEM);
  TString fOutName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  //Checking that files exist
  if(gSystem->AccessPathName(fNamedRICH)){
    cout <<"[ERROR] dRICH run not found\n";
    exit(EXIT_FAILURE);
  }
  if(gSystem->AccessPathName(fNameGEM)){
    cout <<"[ERROR] GEM run not found\n";
    exit(EXIT_FAILURE);
  }
  //Getting the TTrees
  TFile *fdRICH = new TFile(fNamedRICH,"READ");
  if(fdRICH->IsZombie()){
    cout <<"[ERROR] fdRICH is zombie! Something was wrong\n";
    exit(EXIT_FAILURE);
  }
  TTree *t = (TTree*) fdRICH->Get("dRICH");
  t->SetTitle("dRICH and GEM data");
  TFile *fGEM = new TFile(fNameGEM,"READ");
  if(fGEM->IsZombie()){
    cout <<"[ERROR] fGEM is zombie! Something was wrong\n";
    exit(EXIT_FAILURE);
  }
  TTree *tGEM = (TTree*) fGEM->Get("gtr");

  uint nedge;
  int evt, board[MAXDATA], chip[MAXDATA], pol[MAXDATA], time[MAXDATA], marocCh[MAXDATA], pmt[MAXDATA];
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA], ntime[MAXDATA];
  double trigtime;
  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("board",&board);
  t->SetBranchAddress("chip",&chip);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("marocCh",&marocCh);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("radius",&radius);
  t->SetBranchAddress("ntime",&ntime);

  int row, inst, evtGEM;
  float x0, y0, cx0, cy0, x1, y1, cx1, cy1;
  tGEM->SetBranchAddress("row",&row);
  tGEM->SetBranchAddress("inst",&inst);
  tGEM->SetBranchAddress("evt",&evtGEM);
  tGEM->SetBranchAddress("x0",&x0);
  tGEM->SetBranchAddress("y0",&y0);
  tGEM->SetBranchAddress("cx0",&cx0);
  tGEM->SetBranchAddress("cy0",&cy0);
  tGEM->SetBranchAddress("x1",&x1);
  tGEM->SetBranchAddress("y1",&y1);
  tGEM->SetBranchAddress("cx1",&cx1);
  tGEM->SetBranchAddress("cy1",&cy1);

  float gx0=0, gy0=0, gx1=0, gy1=0, gxa=0, gya=0, gxtheta=0, gytheta=0;
  
  TFile *fOut = new TFile(fOutName,"RECREATE"); 
  TTree *tout = new TTree("dRICH","dRICH integrated data");
  
  int tevt, tnedge, tpol[MAXDATA], tboard[MAXDATA], tchip[MAXDATA],ttime[MAXDATA],tpmt[MAXDATA], tchannel[MAXDATA];
  double ttrigtime, tx[MAXDATA],ty[MAXDATA],tr[MAXDATA], tnt[MAXDATA];
  bool dataGEM;
  auto tEvt=tout->Branch("evt",&tevt,"evt/I");
  auto tTrigtime=tout->Branch("trigtime",&ttrigtime,"trigtime/I");
  auto tNedge=tout->Branch("nedge",&tnedge,"nedge/I");
  auto tPol=tout->Branch("pol",&tpol,"pol[nedge]/I");
  auto tTime=tout->Branch("time",&ttime,"time[nedge]/I");
  auto tPmt=tout->Branch("pmt",&tpmt,"pmt[nedge]/I");
  auto tBoard=tout->Branch("board",&tboard,"board[nedge]/I");
  auto tChip=tout->Branch("chip",&tchip,"chip[nedge]/I");
  auto tX=tout->Branch("x",&tx,"x[nedge]/D");
  auto tY=tout->Branch("y",&ty,"y[nedge]/D");
  auto tR=tout->Branch("r",&tr,"r[nedge]/D");
  auto tNT=tout->Branch("nt",&tnt,"nt[nedge]/D");

  TH1D *hX0 = new TH1D("hX0","hX0",200,-100,100);
  TH1D *hY0 = new TH1D("hY0","hY0",200,-100,100);
  TH1D *hX1 = new TH1D("hX1","hX1",200,-100,100);
  TH1D *hY1 = new TH1D("hY1","hY1",200,-100,100);
  vector<int> vGEMentry;

  cout <<"Computing tracking data\n";
  int startCycle = 0;
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    dataGEM=false;
    t->GetEntry(i);
    for(int j = startCycle; j < tGEM->GetEntries(); j++){
      tGEM->GetEntry(j);
      if(evt == evtGEM+1){
        startCycle = max(0,j - 2);
        dataGEM=true;
        vGEMentry.push_back(j);
        break;
      }
      if(evt < evtGEM+1){
        startCycle = max(0,j - 2);
        break; 
      }
    }
    if(dataGEM==false) continue;
    //Copy the dRICH info
    tevt=(int)evt;
    ttrigtime=trigtime;
    tnedge=(int)nedge;
    for(int j = 0; j < nedge; j++){
      tpol[j]=pol[j];
      tboard[j]=board[j];
      tchip[j]=chip[j];
      ttime[j]=time[j];
      tpmt[j]=pmt[j];
      tx[j]=x[j];
      ty[j]=y[j];
      tr[j]=radius[j];
      tnt[j]=ntime[j];
    }
    //Compute GEM info
    float tmpx0=x0;
    float tmpy0=y0;
    float tmpx1=x1;
    float tmpy1=y1;
    GEM_relative(&tmpx0,&tmpy0,&tmpx1,&tmpy1);
    hX0->Fill(tmpx0);
    hY0->Fill(tmpy0);
    hX1->Fill(tmpx1);
    hY1->Fill(tmpy1);
    tout->Fill();
  }
  printEnd();
  run->UpGEMxRunOff=GEM_getBeamlineOffset(hX0);
  run->UpGEMyRunOff=GEM_getBeamlineOffset(hY0);
  run->DnGEMxRunOff=GEM_getBeamlineOffset(hX1);
  run->DnGEMyRunOff=GEM_getBeamlineOffset(hY1);

  auto tGX0=tout->Branch("gx0",&gx0,"gx0/F");
  auto tGY0=tout->Branch("gy0",&gy0,"gy0/F");
  auto tGX1=tout->Branch("gx1",&gx1,"gx1/F");
  auto tGY1=tout->Branch("gy1",&gy1,"gy1/F");
  auto tGXA=tout->Branch("gxa",&gxa,"gxa/F");
  auto tGYA=tout->Branch("gya",&gya,"gya/F");
  auto tGXtheta=tout->Branch("gxtheta",&gxtheta,"gxtheta/F");
  auto tGYtheta=tout->Branch("gytheta",&gytheta,"gytheta/F");

  cout <<"Prototype and tracking data integration\n";
  for(int i = 0; i < tout->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    tout->GetEntry(i);
    tGEM->GetEntry(vGEMentry[i]);
    gx0=x0;
    gy0=y0;
    gx1=x1;
    gy1=y1;
    GEM_relative(&gx0,&gy0,&gx1,&gy1);
    GEM_position(run,&gx0,&gy0,&gx1,&gy1);
    AERO_computing(run,&gxa,&gya,&gxtheta,&gytheta,gx0,gy0,gx1,gy1);
    tGX0->Fill();
    tGY0->Fill();
    tGX1->Fill();
    tGY1->Fill();
    tGXA->Fill();
    tGYA->Fill();
    tGXtheta->Fill();
    tGYtheta->Fill();
  }
  printEnd();

  tout->Write();
  fOut->Close();
  fdRICH->Close();
  fGEM->Close();
}
