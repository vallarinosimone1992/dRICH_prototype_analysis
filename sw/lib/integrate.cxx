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
#include "tracking.h"
#include "getChannel.h"
#include "photoDetPosition.h"
#include "fillMAPS.h"
#include "definition.h"

using namespace std;

map<string,int> map_MAPMT1;
map<string,int> map_MAPMT2;

map<string,int>::iterator it_map_MAPMT1;
map<string,int>::iterator it_map_MAPMT2;


void TTreeIntegration(THeader *runHead){
  int runDRICH=runHead->runNum;
  int runGEM=runHead->runNumGEM; 
  string phDet=runHead->sensor;
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
    cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
    exit(EXIT_FAILURE);
  }
  TString fNamedRICH=Form("%s/processed_data/firstStepData/run_%04d_processed.root",&env_var[0],runDRICH);
  TString fNameGEM=Form("%s/DATA/GEM_DATA/run_%04d_gem.root",&env_var[0],runGEM);
  TString fOutName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&env_var[0],runDRICH);
  cout <<Form("-->> dRICH run: %s\n",&fNamedRICH[0]);
  //Checking that files exist
  if(gSystem->AccessPathName(fNamedRICH)){
    cout <<"[ERROR] dRICH run not found\n";
    exit(EXIT_FAILURE);
  }
  cout <<Form("-->> GEM run: %s\n",&fNameGEM[0]);
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


  cout <<Form("TTrees has %lld & %lld entries\n",t->GetEntries(), tGEM->GetEntries());

  uint nedge;
  int evt, board[MAXDATA], chip[MAXDATA], pol[MAXDATA], time[MAXDATA], marocCh[MAXDATA], pmt[MAXDATA];
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA];
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

  cout <<"dRICH variables setted\n";

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

  cout <<"GEM variables setted\n";

  //int gRow, gInst;
  //float gX0, gY0, gCx0, gCy0, gX1, gY1, gCx1, gCy1, gXup, gYup, gXdown, gYdown, gXaero, gYaero, angleX, angleY, gRelXup, gRelYup, gRelXdown, gRelYdown; 
  float gx0=0, gy0=0, gx1=0, gy1=0, gxa=0, gya=0, gxtheta=0, gytheta=0;
  
  TFile *fOut = new TFile(fOutName,"RECREATE"); 
  TTree *tout = new TTree("dRICH","dRICH integrated data");
  
  int tevt, tnedge, tpol[MAXDATA], tboard[MAXDATA], tchip[MAXDATA],ttime[MAXDATA],tpmt[MAXDATA], tchannel[MAXDATA];
  double ttrigtime, tx[MAXDATA],ty[MAXDATA],tr[MAXDATA];
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

  TH1D *hX0 = new TH1D("hX0","hX0",200,-100,100);
  TH1D *hY0 = new TH1D("hY0","hY0",200,-100,100);
  TH1D *hX1 = new TH1D("hX1","hX1",200,-100,100);
  TH1D *hY1 = new TH1D("hY1","hY1",200,-100,100);

  vector<int> vGEMentry;


  int startCycle = 0;
  for(int i = 0; i < t->GetEntries(); i++){
    if((i)%(t->GetEntries()/10)==0)cout <<Form("\rSynchronizing the dRich and GEM data: %lld%% completed    ",(1+(i/(t->GetEntries()/10)))*10) <<flush;
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
    //cout <<tnedge <<" " <<nedge <<endl;
    for(int j = 0; j < nedge; j++){
      tpol[j]=pol[j];
      tboard[j]=board[j];
      tchip[j]=chip[j];
      ttime[j]=time[j];
      tpmt[j]=pmt[j];
      tx[j]=x[j];
      ty[j]=y[j];
      tr[j]=radius[j];
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
    /*
    gx0=x0;
    gy0=y0;
    gx1=x1;
    gy1=y1;
    GEM_relative(&gx0,&gy0,&gx1,&gy1);
    hX0->Fill(gx0);
    hY0->Fill(gy0);
    hX1->Fill(gx1);
    hY1->Fill(gy1);
    tout->Fill();*/
  }
  cout <<endl;
  //tout->Write();

  runHead->UpGEMxRunOff=GEM_getBeamlineOffset(hX0);
  runHead->UpGEMyRunOff=GEM_getBeamlineOffset(hY0);
  runHead->DnGEMxRunOff=GEM_getBeamlineOffset(hX1);
  runHead->DnGEMyRunOff=GEM_getBeamlineOffset(hY1);

  cout <<"GEM z: " <<runHead->UpGEMz <<" " <<runHead->DnGEMz <<endl;
  cout <<"Aero z: " <<runHead->zAerogel <<endl;
  cout <<"GEM entry size " <<vGEMentry.size() <<" " <<tout->GetEntries() <<endl;
  auto tGX0=tout->Branch("gx0",&gx0,"gx0/F");
  auto tGY0=tout->Branch("gy0",&gy0,"gy0/F");
  auto tGX1=tout->Branch("gx1",&gx1,"gx1/F");
  auto tGY1=tout->Branch("gy1",&gy1,"gy1/F");
  auto tGXA=tout->Branch("gxa",&gxa,"gxa/F");
  auto tGYA=tout->Branch("gya",&gya,"gya/F");
  auto tGXtheta=tout->Branch("gxtheta",&gxtheta,"gxtheta/F");
  auto tGYtheta=tout->Branch("gytheta",&gytheta,"gytheta/F");

  for(int i = 0; i < tout->GetEntries(); i++){
    if((i)%(tout->GetEntries()/10)==0)cout <<Form("\rComputing tracking info: %lld%% completed    ",(1+(i/(tout->GetEntries()/10)))*10) <<flush;
    tout->GetEntry(i);
    tGEM->GetEntry(vGEMentry[i]);
    gx0=x0;
    gy0=y0;
    gx1=x1;
    gy1=y1;
    GEM_relative(&gx0,&gy0,&gx1,&gy1);
    //cout <<"Relative: " <<gx0 <<" " <<gy0 <<" " <<gx1 <<" " <<gy1 <<endl;
    GEM_position(runHead,&gx0,&gy0,&gx1,&gy1);
    AERO_computing(runHead,&gxa,&gya,&gxtheta,&gytheta,gx0,gy0,gx1,gy1);
    //cout <<"Check here: "<<gxa <<" " <<gya <<" " <<gxtheta <<" " <<gytheta <<endl;
    tGX0->Fill();
    tGY0->Fill();
    tGX1->Fill();
    tGY1->Fill();
    tGXA->Fill();
    tGYA->Fill();
    tGXtheta->Fill();
    tGYtheta->Fill();
  }
  cout <<endl;

  tout->Write();
  fOut->Close();
  fdRICH->Close();
  fGEM->Close();
  cout <<"End of the function\n";
}


//Add header
/* auto tHeader=tout->Branch("header",THeader,&runHead);
   (void) tHeader; 
   auto tRunNum=tout->Branch("runNum",&runHead->runNum,"runNum/I");
   (void) tRunNum;
   auto tEnergyGeV=tout->Branch("energyGeV",&runHead->energyGeV,"energyGeV/I");
   (void) tEnergyGeV;
   auto tExpEvents=tout->Branch("expEvents",&runHead->expEvents,"expEvents/I");
   (void) tExpEvents;
   auto tPowerHV=tout->Branch("powerHV",&runHead->powerHV,"powerHV/I");
   (void) tPowerHV;
   auto tRunGEM=tout->Branch("runNumGEM",&runHead->runNumGEM,"runNumGEM/I");
   (void) tRunGEM;
   auto tPedestalGEM=tout->Branch("pedestalGEM",&runHead->pedestalGEM,"pedestalGEM/I");
   (void) tPedestalGEM;

   auto tDay=tout->Branch("day",&runHead->day,"day/C");
   (void) tDay;
   */




/*
   void TTreeIntegration(THeader *runHead){
   int runDRICH=runHead->runNum;
   int runGEM=runHead->runNumGEM; 
   string phDet=runHead->sensor;
   const char  *tmp = getenv("DRICH_SUITE");
   getMapMAPMT(&map_MAPMT1,&map_MAPMT2);
   string env_var(tmp ? tmp : "");
   if(env_var.empty()){
   cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
   exit(EXIT_FAILURE);
   }
   TString fNamedRICH=Form("%s/DATA/dRICH_DATA/run_%04d.root",&env_var[0],runDRICH);
   TString fNameGEM=Form("%s/DATA/GEM_DATA/run_%04d_gem.root",&env_var[0],runGEM);
   TString fOutName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_processed.root",&env_var[0],runDRICH);
   cout <<Form("-->> dRICH run: %s\n",&fNamedRICH[0]);
   if(gSystem->AccessPathName(fNamedRICH)){
   cout <<"[ERROR] dRICH run not found\n";
   exit(EXIT_FAILURE);
   }
   cout <<Form("-->> GEM run: %s\n",&fNameGEM[0]);
   if(gSystem->AccessPathName(fNameGEM)){
   cout <<"[ERROR] GEM run not found\n";
   exit(EXIT_FAILURE);
   }
   TFile *fdRICH = new TFile(fNamedRICH,"READ");
   if(fdRICH->IsZombie()){
   cout <<"[ERROR] fdRICH is zombie! Something was wrong\n";
   exit(EXIT_FAILURE);
   }
   TFile *fGEM = new TFile(fNameGEM,"READ");
   if(fGEM->IsZombie()){
   cout <<"[ERROR] fGEM is zombie! Something was wrong\n";
   exit(EXIT_FAILURE);
   }

   TTree *t = (TTree*) fdRICH->Get("data");
   TTree *tGEM = (TTree*) fGEM->Get("gtr");

   cout <<Form("TTrees has %lld & %lld entries\n",t->GetEntries(), tGEM->GetEntries());

   TFile *fOut = new TFile(fOutName,"RECREATE"); 
   TTree *tout = new TTree("dRICH","dRICH integrated data");

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

int evt, tSlot[MAXDATA], tFiber[MAXDATA], tCh[MAXDATA], tPol[MAXDATA],tTime[MAXDATA];
uint tNedge;
double tTrigTime;
double x[MAXDATA], y[MAXDATA], radius[MAXDATA];
int pmt[MAXDATA], channel[MAXDATA];
cout <<"GEM variables setted\n";

//Add header
auto tHeader=tout->Branch("header",THeader,&runHead);
(void) tHeader; 
auto tRunNum=tout->Branch("runNum",&runHead->runNum,"runNum/I");
(void) tRunNum;
auto tEnergyGeV=tout->Branch("energyGeV",&runHead->energyGeV,"energyGeV/I");
(void) tEnergyGeV;
auto tExpEvents=tout->Branch("expEvents",&runHead->expEvents,"expEvents/I");
(void) tExpEvents;
auto tPowerHV=tout->Branch("powerHV",&runHead->powerHV,"powerHV/I");
(void) tPowerHV;
auto tRunGEM=tout->Branch("runNumGEM",&runHead->runNumGEM,"runNumGEM/I");
(void) tRunGEM;
auto tPedestalGEM=tout->Branch("pedestalGEM",&runHead->pedestalGEM,"pedestalGEM/I");
(void) tPedestalGEM;

auto tDay=tout->Branch("day",&runHead->day,"day/C");
(void) tDay;

auto evento=tout->Branch("evt",&evt,"evt/I");
(void)evento;
auto ttrigtime=tout->Branch("trigtime",&tTrigTime,"trigtime/D");
(void)ttrigtime;
auto tnedge=tout->Branch("nedge",&tNedge,"nedge/i");
(void)tnedge;
auto tslot=tout->Branch("slot",&tSlot[0],"slot[nedge]/I");
(void)tslot;
auto tfiber=tout->Branch("fiber",&tFiber[0],"fiber[nedge]/I");
(void)tfiber;
auto tch=tout->Branch("ch",&tCh[0],"ch[nedge]/I");
(void)tch;
auto tpol=tout->Branch("pol",&tPol[0],"pol[nedge]/I");
(void)tpol;
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

int gRow, gInst;
float gX0, gY0, gCx0, gCy0, gX1, gY1, gCx1, gCy1, gXup, gYup, gXdown, gYdown, gXaero, gYaero, angleX, angleY, gRelXup, gRelYup, gRelXdown, gRelYdown;


int startCycle = 0;
for(int i = 0; i < t->GetEntries(); i++){
  if((i)%(t->GetEntries()/10)==0)cout <<Form("\rSynchronizing the dRich and GEM data: %lld%% completed    ",(1+(i/(t->GetEntries()/10)))*10) <<flush;
  bool flagEvt = false;
  t->GetEntry(i);
  for(int j = startCycle; j < tGEM->GetEntries(); j++){
    tGEM->GetEntry(j);
    if(evtDRICH == (unsigned int) evtGEM+1){
      evt = evtDRICH;
      startCycle = max(0,j - 2);
      flagEvt = true;
      break;
    }
  } 
  if(flagEvt == false) continue;
  //Copy the dRich data
  tTrigTime=trigtime;
  tNedge=nedge;
  for(uint j = 0; j < nedge; j++){
    tSlot[j]=slot[j];
    tFiber[j]=fiber[j];
    tCh[j]=chM[j];
    tPol[j]=pol[j];
    tTime[j]=time[j];
    //if(pol[j]==0) hTime->Fill(time[j]);
    //Producing the position, radius, pmt and channel info
    if(phDet.compare("MAPMT")==0){
      if(fiber[j] == 4 || fiber[j] == 5 || fiber[j] == 6 || fiber[j] == 7 ){ 
        channel[j]= map_MAPMT1.at(getMAPMT_ch(fiber[j],chM[j]));
      }else if(fiber[j] == 8 || fiber[j] == 9 || fiber[j] == 10 || fiber[j] == 11 ){
        channel[j]= map_MAPMT2.at(getMAPMT_ch(fiber[j],chM[j]));
      }else{continue;}
    }
    pmt[j]=FiberToPlace(fiber[j]);
    //if(phDet.compare("MPPC")==0)MPPCposition(channel[j],pmt[j],&x[j],&y[j]);
    if(phDet.compare("MAPMT")==0)MAPMTposition(channel[j],pmt[j],&x[j],&y[j]);
    radius[j]=sqrt(pow(x[j],2)+pow(y[j],2));
  }
  //Copy the GEM info
  gRow=row;
  gInst=inst;
  gX0=x0;
  gY0=y0;
  gCx0=cx0;
  gCy0=cy0;
  gX1=x1;
  gY1=y1;
  gCx1=cx1;
  gCy1=cy1;
  tout->Fill();
}
cout <<endl;
tout->Write();
fOut->Close();
cout <<"End of the function\n";
}
*/
