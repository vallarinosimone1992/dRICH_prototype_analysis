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
#include <TStyle.h>
#include <TCanvas.h>

#include "integrate.h"
#include "utility.h"
#include "tracking.h"
#include "getChannel.h"
#include "photoDetPosition.h"
#include "fillMAPS.h"
#include "definition.h"
#include "computing.h"
#include "writeHeaderText.h"

using namespace std;

void prepareMergedRun(THeader *run){
  TString fNamedRICH=Form("%s/DATA/MERGER_DATA/merged_run_%04d.root",run->suite.c_str(),run->runNum);
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
  TTree *t = (TTree*) fdRICH->Get("kTarget_dRICH");
  t->SetTitle("dRICH and GEM data");

  int nedge;
  int evt, board[MAXDATA], chip[MAXDATA], pol[MAXDATA], slot[MAXDATA], fiber[MAXDATA], ch[MAXDATA], time[MAXDATA], anode[MAXDATA], pmt[MAXDATA], otime[MAXDATA];
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA], ntime[MAXDATA];
  double trigtime, x474time, x519time, x537time;
  bool trigSig[MAXDATA];

  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("slot",&slot);
  t->SetBranchAddress("fiber",&fiber);
  t->SetBranchAddress("ch",&ch);
  t->SetBranchAddress("board",&board);
  t->SetBranchAddress("chip",&chip);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("origTime",&otime);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("anode",&anode);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("x474time",&x474time);
  t->SetBranchAddress("x519time",&x519time);
  t->SetBranchAddress("x537time",&x537time);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("radius",&radius);
  t->SetBranchAddress("nt",&ntime);
  t->SetBranchAddress("trigSig",&trigSig);
  
  int row, inst, evtGEM;
  float x0, y0, cx0, cy0, x1, y1, cx1, cy1;
  t->SetBranchAddress("row",&row);
  t->SetBranchAddress("inst",&inst);
  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("x0",&x0);
  t->SetBranchAddress("y0",&y0);
  t->SetBranchAddress("cx0",&cx0);
  t->SetBranchAddress("cy0",&cy0);
  t->SetBranchAddress("x1",&x1);
  t->SetBranchAddress("y1",&y1);
  t->SetBranchAddress("cx1",&cx1);
  t->SetBranchAddress("cy1",&cy1);

  float gx0=0, gy0=0, gx1=0, gy1=0, gxa=0, gya=0, gxtheta=0, gytheta=0;
  
  TFile *fOut = new TFile(fOutName,"RECREATE"); 
  TTree *tout = new TTree("dRICH","dRICH integrated data");
  
  int tevt, tnedge, tpol[MAXDATA], tslot[MAXDATA], tfiber[MAXDATA], tch[MAXDATA], tboard[MAXDATA], tchip[MAXDATA],ttime[MAXDATA],tpmt[MAXDATA], tchannel[MAXDATA], totime[MAXDATA];
  double ttrigtime, tx474time, tx519time, tx537time, tx[MAXDATA],ty[MAXDATA],trMM[MAXDATA], tr[MAXDATA], tnt[MAXDATA];
  bool dataGEM, ttrigSig[MAXDATA];
  float tx0, tx1, ty0, ty1;

  auto tEvt=tout->Branch("evt",&tevt,"evt/I");
  auto tTrigtime=tout->Branch("trigtime",&ttrigtime,"trigtime/D");
  auto tX474time=tout->Branch("x474time",&tx474time,"x474time/D");
  auto tX519time=tout->Branch("x519time",&tx519time,"x519time/D");
  auto tX537time=tout->Branch("x537time",&tx537time,"x537time/D");
  auto tNedge=tout->Branch("nedge",&tnedge,"nedge/I");
  auto tPol=tout->Branch("pol",&tpol,"pol[nedge]/I");
  auto tSlot=tout->Branch("slot",&tslot,"slot[nedge]/I");
  auto tFiber=tout->Branch("fiber",&tfiber,"fiber[nedge]/I");
  auto tCh=tout->Branch("ch",&tch,"ch[nedge]/I");
  auto tTime=tout->Branch("time",&ttime,"time[nedge]/I");
  auto toTime=tout->Branch("otime",&totime,"original_time[nedge]/I");
  auto tPmt=tout->Branch("pmt",&tpmt,"pmt[nedge]/I");
  auto tBoard=tout->Branch("board",&tboard,"board[nedge]/I");
  auto tChip=tout->Branch("chip",&tchip,"chip[nedge]/I");
  auto tX=tout->Branch("x",&tx,"x[nedge]/D");
  auto tY=tout->Branch("y",&ty,"y[nedge]/D");
  auto tR=tout->Branch("r",&tr,"r[nedge]/D");
  auto tRmm=tout->Branch("rmm",&trMM,"rmm[nedge]/D");
  auto tNT=tout->Branch("nt",&tnt,"nt[nedge]/D");
  auto tTrigSig=tout->Branch("trigSig",&ttrigSig,"trigSig[nedge]/O");

  
  auto tGX0=tout->Branch("gx0",&gx0,"gx0/F");
  auto tGY0=tout->Branch("gy0",&gy0,"gy0/F");
  auto tGX1=tout->Branch("gx1",&gx1,"gx1/F");
  auto tGY1=tout->Branch("gy1",&gy1,"gy1/F");
  auto tGXA=tout->Branch("gxa",&gxa,"gxa/F");
  auto tGYA=tout->Branch("gya",&gya,"gya/F");
  auto tGXtheta=tout->Branch("gxtheta",&gxtheta,"gxtheta/F");
  auto tGYtheta=tout->Branch("gytheta",&gytheta,"gytheta/F");
  
  
  cout <<"Taking one of the run header\n";
  cout <<"The run number is " <<run->runNumGEM <<endl;

  THeader firstRunHeader;
  //readHeaders(run->runNumGEM,&firstRunHeader);
  readHeaders(run->pedestalGEM,&firstRunHeader);
  readHeaderShort(&firstRunHeader);

  run->UpGEMxRunOff=firstRunHeader.UpGEMxRunOff;
  run->UpGEMyRunOff=firstRunHeader.UpGEMyRunOff;
  run->DnGEMxRunOff=firstRunHeader.DnGEMxRunOff;
  run->DnGEMyRunOff=firstRunHeader.DnGEMyRunOff;


  cout <<"No tracking info. I'm filling the tracking data with null values\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    dataGEM=true;
    t->GetEntry(i);
    //Copy the dRICH info
    tevt=(int)evt;
    ttrigtime=trigtime;
    tx474time=x474time;
    tx519time=x519time;
    tx537time=x537time;
    tnedge=nedge;
    for(int j = 0; j < nedge; j++){
      tpol[j]=pol[j];
      tslot[j]=slot[j];
      tfiber[j]=fiber[j];
      tch[j]=anode[j];
      tboard[j]=board[j];
      tchip[j]=chip[j];
      ttime[j]=time[j];
      totime[j]=otime[j];
      tpmt[j]=pmt[j];
      tx[j]=x[j];
      ty[j]=y[j];
      double inPath=0;
      double zMir=0;
      if(radius[j] > run->geoCut){
        inPath=run->firstPath;
        zMir=run->firstMirrorPosition;
      }else{
        inPath=run->secondPath;
        zMir=run->secondMirrorPosition;
      }
      trMM[j]=radius[j];
      tr[j]=mmTomRad(radius[j],inPath,zMir);
      tnt[j]=ntime[j];
      ttrigSig[j]=trigSig[j];
    }
    if(SWAP_UPSTREAM_DOWNSTREAM_GEM){
      gx0=x1;
      gy0=y1;
      gx1=x0;
      gy1=y0;
    }else{
      gx0=x0;
      gy0=y0;
      gx1=x1;
      gy1=y1;
    }
    GEM_relative(&gx0,&gy0,&gx1,&gy1);
    GEM_position(run,&gx0,&gy0,&gx1,&gy1);
    AERO_computing(run,&gxa,&gya,&gxtheta,&gytheta,gx0,gy0,gx1,gy1);

    tout->Fill();
  }
  printEnd();

  tout->Write();
  fOut->Close();
  fdRICH->Close();
}


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
  int evt, board[MAXDATA], chip[MAXDATA], pol[MAXDATA], slot[MAXDATA], fiber[MAXDATA], ch[MAXDATA], time[MAXDATA], anode[MAXDATA], pmt[MAXDATA], otime[MAXDATA];
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA], ntime[MAXDATA];
  double trigtime, x474time, x519time, x537time;
  bool trigSig[MAXDATA];

  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("x474time",&x474time);
  t->SetBranchAddress("x519time",&x519time);
  t->SetBranchAddress("x537time",&x537time);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("slot",&slot);
  t->SetBranchAddress("fiber",&fiber);
  t->SetBranchAddress("ch",&ch);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("origTime",&otime);
  t->SetBranchAddress("board",&board);
  t->SetBranchAddress("chip",&chip);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("anode",&anode);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("radius",&radius);
  t->SetBranchAddress("nt",&ntime);
  t->SetBranchAddress("trigSig",&trigSig);


  float gx0=0, gy0=0, gx1=0, gy1=0, gxa=0, gya=0, gxtheta=0, gytheta=0;

  TFile *fOut = new TFile(fOutName,"RECREATE"); 
  TTree *tout = new TTree("dRICH","dRICH integrated data");

  int tevt, tnedge, tpol[MAXDATA], tslot[MAXDATA], tfiber[MAXDATA], tch[MAXDATA], tboard[MAXDATA], tchip[MAXDATA],ttime[MAXDATA],tpmt[MAXDATA], tchannel[MAXDATA], totime[MAXDATA];
  double ttrigtime, tx474time, tx519time, tx537time, tx[MAXDATA],ty[MAXDATA],trMM[MAXDATA], tr[MAXDATA], tnt[MAXDATA];
  bool dataGEM, ttrigSig[MAXDATA];

  auto tEvt=tout->Branch("evt",&tevt,"evt/I");
  auto tTrigtime=tout->Branch("trigtime",&ttrigtime,"trigtime/D");
  auto tX474time=tout->Branch("x474time",&tx474time,"x474time/D");
  auto tX519time=tout->Branch("x519time",&tx519time,"x519time/D");
  auto tX537time=tout->Branch("x537time",&tx537time,"x537time/D");
  auto tNedge=tout->Branch("nedge",&tnedge,"nedge/I");
  auto tPol=tout->Branch("pol",&tpol,"pol[nedge]/I");
  auto tSlot=tout->Branch("slot",&tslot,"slot[nedge]/I");
  auto tFiber=tout->Branch("fiber",&tfiber,"fiber[nedge]/I");
  auto tCh=tout->Branch("ch",&tch,"ch[nedge]/I");
  auto tTime=tout->Branch("time",&ttime,"time[nedge]/I");
  auto toTime=tout->Branch("otime",&totime,"original_time[nedge]/I");
  auto tPmt=tout->Branch("pmt",&tpmt,"pmt[nedge]/I");
  auto tBoard=tout->Branch("board",&tboard,"board[nedge]/I");
  auto tChip=tout->Branch("chip",&tchip,"chip[nedge]/I");
  auto tX=tout->Branch("x",&tx,"x[nedge]/D");
  auto tY=tout->Branch("y",&ty,"y[nedge]/D");
  auto tR=tout->Branch("r",&tr,"r[nedge]/D");
  auto tRmm=tout->Branch("rmm",&trMM,"rmm[nedge]/D");
  auto tNT=tout->Branch("nt",&tnt,"nt[nedge]/D");
  auto tTrigSig=tout->Branch("trigSig",&ttrigSig,"trigSig[nedge]/O");

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
    tx474time=x474time;
    tx519time=x519time;
    tx537time=x537time;
    tnedge=(int)nedge;
    for(int j = 0; j < nedge; j++){
      tpol[j]=pol[j];
      tslot[j]=slot[j];
      tfiber[j]=fiber[j];
      tch[j]=anode[j];
      tboard[j]=board[j];
      tchip[j]=chip[j];
      ttime[j]=time[j];
      totime[j]=otime[j];
      tpmt[j]=pmt[j];
      tx[j]=x[j];
      ty[j]=y[j];
      double inPath=0;
      double zMir=0;
      if(radius[j] > run->geoCut){
        inPath=run->firstPath;
        zMir=run->firstMirrorPosition;
      }else{
        inPath=run->secondPath;
        zMir=run->secondMirrorPosition;
      }
      trMM[j]=radius[j];
      tr[j]=mmTomRad(radius[j],inPath,zMir);
      tnt[j]=ntime[j];
      ttrigSig[j]=trigSig[j];
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
  int evt, board[MAXDATA], chip[MAXDATA], pol[MAXDATA], slot[MAXDATA], fiber[MAXDATA], ch[MAXDATA], time[MAXDATA], anode[MAXDATA], pmt[MAXDATA], otime[MAXDATA];
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA], ntime[MAXDATA];
  double trigtime, x474time, x519time, x537time;
  bool trigSig[MAXDATA];

  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("x474time",&x474time);
  t->SetBranchAddress("x519time",&x519time);
  t->SetBranchAddress("x537time",&x537time);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("slot",&slot);
  t->SetBranchAddress("fiber",&fiber);
  t->SetBranchAddress("ch",&ch);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("origTime",&otime);
  t->SetBranchAddress("board",&board);
  t->SetBranchAddress("chip",&chip);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("anode",&anode);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("radius",&radius);
  t->SetBranchAddress("nt",&ntime);
  t->SetBranchAddress("trigSig",&trigSig);

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

  int tevt, tnedge, tpol[MAXDATA], tslot[MAXDATA], tfiber[MAXDATA], tch[MAXDATA], tboard[MAXDATA], tchip[MAXDATA],ttime[MAXDATA],tpmt[MAXDATA], tchannel[MAXDATA], totime[MAXDATA];
  double ttrigtime, tx474time, tx519time, tx537time, tx[MAXDATA],ty[MAXDATA],trMM[MAXDATA], tr[MAXDATA], tnt[MAXDATA];
  bool dataGEM, ttrigSig[MAXDATA];
  auto tEvt=tout->Branch("evt",&tevt,"evt/I");
  auto tTrigtime=tout->Branch("trigtime",&ttrigtime,"trigtime/D");
  auto tX474time=tout->Branch("x474time",&tx474time,"x474time/D");
  auto tX519time=tout->Branch("x519time",&tx519time,"x519time/D");
  auto tX537time=tout->Branch("x537time",&tx537time,"x537time/D");
  auto tNedge=tout->Branch("nedge",&tnedge,"nedge/I");
  auto tPol=tout->Branch("pol",&tpol,"pol[nedge]/I"); 
  auto tSlot=tout->Branch("slot",&tslot,"slot[nedge]/I");
  auto tFiber=tout->Branch("fiber",&tfiber,"fiber[nedge]/I");
  auto tCh=tout->Branch("ch",&tch,"ch[nedge]/I");
  auto tTime=tout->Branch("time",&ttime,"time[nedge]/I");
  auto toTime=tout->Branch("otime",&totime,"original_time[nedge]/I");
  auto tPmt=tout->Branch("pmt",&tpmt,"pmt[nedge]/I");
  auto tBoard=tout->Branch("board",&tboard,"board[nedge]/I");
  auto tChip=tout->Branch("chip",&tchip,"chip[nedge]/I");
  auto tX=tout->Branch("x",&tx,"x[nedge]/D");
  auto tY=tout->Branch("y",&ty,"y[nedge]/D");
  auto tR=tout->Branch("r",&tr,"r[nedge]/D");
  auto tRmm=tout->Branch("rmm",&trMM,"rmm[nedge]/D");
  auto tNT=tout->Branch("nt",&tnt,"nt[nedge]/D");
  auto tTrigSig=tout->Branch("trigSig",&ttrigSig,"trigSig[nedge]/O");

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
      if(evt < evtGEM+2){
        startCycle = max(0,j - 2);
        break; 
      }
    }
    if(dataGEM==false) continue;
    //Copy the dRICH info
    tevt=(int)evt;
    ttrigtime=trigtime;
    tx474time=x474time;
    tx519time=x519time;
    tx537time=x537time;
    tnedge=(int)nedge;
    for(int j = 0; j < nedge; j++){
      tpol[j]=pol[j];
      tslot[j]=slot[j];
      tfiber[j]=fiber[j];
      tch[j]=anode[j];
      tboard[j]=board[j];
      tchip[j]=chip[j];
      ttime[j]=time[j];
      totime[j]=otime[j];
      tpmt[j]=pmt[j];
      tx[j]=x[j];
      ty[j]=y[j];
      double inPath=0;
      double zMir=0;
      if(radius[j] > run->geoCut){
        inPath=run->firstPath;
        zMir=run->firstMirrorPosition;
        //zMir=run->secondMirrorPosition; //SWAP TO CHECK THE GEM
      }else{
        inPath=run->secondPath;
        //zMir=run->firstMirrorPosition;//SWAP TO CHECK THE GEM
        zMir=run->secondMirrorPosition;
      }
      trMM[j]=radius[j];
      tr[j]=mmTomRad(radius[j],inPath,zMir);
      tnt[j]=ntime[j];
      ttrigSig[j]=trigSig[j];
    }
    //Compute GEM info
    float tmpx0=x0;
    float tmpy0=y0;
    float tmpx1=x1;
    float tmpy1=y1;
    if(SWAP_UPSTREAM_DOWNSTREAM_GEM){
      tmpx0=x1;
      tmpy0=y1;
      tmpx1=x0;
      tmpy1=y0;
    }


    GEM_relative(&tmpx0,&tmpy0,&tmpx1,&tmpy1);
    hX0->Fill(tmpx0);
    hY0->Fill(tmpy0);
    hX1->Fill(tmpx1);
    hY1->Fill(tmpy1);
    tout->Fill();
  }
  printEnd();
  run->UpGEMxRunOff=GEM_getBeamlineOffset(hX0,0);
  run->UpGEMyRunOff=GEM_getBeamlineOffset(hY0,0);
  run->DnGEMxRunOff=GEM_getBeamlineOffset(hX1,1);
  run->DnGEMyRunOff=GEM_getBeamlineOffset(hY1,0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TCanvas *cGEM = new TCanvas("cGEM","cGEM",1600,900);
  cGEM->Divide(2,2);
  cGEM->cd(1);
  hX0->Draw();
  cGEM->cd(2);
  hY0->Draw();
  cGEM->cd(3);
  hX1->Draw();
  cGEM->cd(4);
  hY1->Draw();
  cGEM->Update();
  //cGEM->Print("cGEM.pdf");
  //cGEM->Print("cGEM.root");
  cGEM->Close();
  /*run->UpGEMxRunOff=0;
    run->UpGEMyRunOff=0;
    run->DnGEMxRunOff=0;
    run->DnGEMyRunOff=0;*/

  auto tGX0=tout->Branch("gx0",&gx0,"gx0/F");
  auto tGY0=tout->Branch("gy0",&gy0,"gy0/F");
  auto tGX1=tout->Branch("gx1",&gx1,"gx1/F");
  auto tGY1=tout->Branch("gy1",&gy1,"gy1/F");
  auto tGXA=tout->Branch("gxa",&gxa,"gxa/F");
  auto tGYA=tout->Branch("gya",&gya,"gya/F");
  auto tGXtheta=tout->Branch("gxtheta",&gxtheta,"gxtheta/F");
  auto tGYtheta=tout->Branch("gytheta",&gytheta,"gytheta/F");

  ///////// ATTENTION! I SWAP THE GEM TO CHECK IT.
  cout <<"Prototype and tracking data integration\n";
  for(int i = 0; i < tout->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    tout->GetEntry(i);
    tGEM->GetEntry(vGEMentry[i]);
    ///////// ATTENTION! I SWAP THE GEM TO CHECK IT.
    if(SWAP_UPSTREAM_DOWNSTREAM_GEM){
      gx0=x1;
      gy0=y1;
      gx1=x0;
      gy1=y0;
    }else{
      gx0=x0;
      gy0=y0;
      gx1=x1;
      gy1=y1;
    }
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
    cout <<Form("%lf %lf %lf %lf\n",gx0,gy0,gx1,gy1);
  }
  printEnd();

  tout->Write();
  fOut->Close();
  fdRICH->Close();
  fGEM->Close();
}
