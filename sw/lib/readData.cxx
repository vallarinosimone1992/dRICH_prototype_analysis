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
    cout <<"[ERROR] Inconsistency between the selected data type and logbook [SENSOR]\n";
    exit(EXIT_FAILURE);
  }
  //Take the MAPMT maps
  map<string,int> map_MAPMT1;
  map<string,int> map_MAPMT2;
  map<string,int>::iterator it_map_MAPMT1;
  map<string,int>::iterator it_map_MAPMT2;
  getMapMAPMT(&map_MAPMT1,&map_MAPMT2);

  //Find DRICH_SUITE environment variable
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
    cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
    exit(EXIT_FAILURE);
  }

  TString fNamedRICH=Form("%s/DATA/dRICH_DATA/run_%04d.root",&env_var[0],runDRICH);
  TString fOutName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_processed.root",&env_var[0],runDRICH);
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

  int evt, tPol[MAXDATA],tTime[MAXDATA];
  uint tNedge;
  double tTrigTime;
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA];
  int pmt[MAXDATA], channel[MAXDATA];


  auto evento=tout->Branch("evt",&evt,"evt/I");
  (void)evento;
  auto ttrigtime=tout->Branch("trigtime",&tTrigTime,"trigtime/D");
  (void)ttrigtime;
  auto tnedge=tout->Branch("nedge",&tNedge,"nedge/i");
  (void)tnedge;
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

  for(int i = 0; i < t->GetEntries(); i++){
    if((i)%(t->GetEntries()/10)==0)cout <<Form("\rSynchronizing the dRich and GEM data: %lld%% completed    ",(1+(i/(t->GetEntries()/10)))*10) <<flush;
    t->GetEntry(i);
    tTrigTime=trigtime;
    tNedge=nedge;
    bool wrongEvent=false;
    for(uint j = 0; j < nedge; j++){
      tPol[j]=pol[j];
      tTime[j]=time[j];
      if(fiber[j] == 4 || fiber[j] == 5 || fiber[j] == 6 || fiber[j] == 7 ){
        channel[j]= map_MAPMT1.at(getMAPMT_ch(fiber[j],chM[j]));
      }else if(fiber[j] == 8 || fiber[j] == 9 || fiber[j] == 10 || fiber[j] == 11 ){
        channel[j]= map_MAPMT2.at(getMAPMT_ch(fiber[j],chM[j]));
      }else{
        wrongEvent=true;
        cout <<"[ERROR] Event with a wrong \"fiber\" value; Analysis will continue without this event\n";
        break;
      }
      pmt[j]=FiberToPhDet(fiber[j],&(runHead->fiberRef)[0]);
      MAPMTposition(channel[j],pmt[j],&x[j],&y[j],&radius[j]);
    }
    if(wrongEvent==true)continue;
    tout->Fill();
  }
  cout <<endl;
  tout->Write();
  fOut->Close();
  cout <<"MAPMTs data wrote\n";

}


// ####################################################### //

void getMPPC(THeader *run){
  int runDRICH=runHead->runNum;
  string phDet=runHead->sensor;
  if(phDet!="MPPC"){
    cout <<"[ERROR] Inconsistency between the selected data type and logbook [SENSOR]\n";
    exit(EXIT_FAILURE);
  }
  //Take the MAPMT maps
  map<string,int> map_MAPMT1;
  map<string,int> map_MAPMT2;
  map<string,int>::iterator it_map_MAPMT1;
  map<string,int>::iterator it_map_MAPMT2;
  getMapMAPMT(&map_MAPMT1,&map_MAPMT2);

  //Find DRICH_SUITE environment variable
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
    cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
    exit(EXIT_FAILURE);
  }

  TString fNamedRICH=Form("%s/DATA/dRICH_DATA/run_%04d.root",&env_var[0],runDRICH);
  TString fOutName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_processed.root",&env_var[0],runDRICH);
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

  int evt, tPol[MAXDATA],tTime[MAXDATA];
  uint tNedge;
  double tTrigTime;
  double x[MAXDATA], y[MAXDATA], radius[MAXDATA];
  int pmt[MAXDATA], channel[MAXDATA];


  auto evento=tout->Branch("evt",&evt,"evt/I");
  (void)evento;
  auto ttrigtime=tout->Branch("trigtime",&tTrigTime,"trigtime/D");
  (void)ttrigtime;
  auto tnedge=tout->Branch("nedge",&tNedge,"nedge/i");
  (void)tnedge;
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

  for(int i = 0; i < t->GetEntries(); i++){
    if((i)%(t->GetEntries()/10)==0)cout <<Form("\rSynchronizing the dRich and GEM data: %lld%% completed    ",(1+(i/(t->GetEntries()/10)))*10) <<flush;
    t->GetEntry(i);
    tTrigTime=trigtime;
    tNedge=nedge;
    bool wrongEvent=false;
    for(uint j = 0; j < nedge; j++){
      tPol[j]=pol[j];
      tTime[j]=time[j];
      if(fiber[j] == 4 || fiber[j] == 5 || fiber[j] == 6 || fiber[j] == 7 ){
        channel[j]= map_MAPMT1.at(getMAPMT_ch(fiber[j],chM[j]));
      }else if(fiber[j] == 8 || fiber[j] == 9 || fiber[j] == 10 || fiber[j] == 11 ){
        channel[j]= map_MAPMT2.at(getMAPMT_ch(fiber[j],chM[j]));
      }else{
        wrongEvent=true;
        cout <<"[ERROR] Event with a wrong \"fiber\" value; Analysis will continue without this event\n";
        break;
      }
      pmt[j]=FiberToMAPMT(fiber[j]);
      MAPMTposition(channel[j],pmt[j],&x[j],&y[j],&radius[j]);
    }
    if(wrongEvent==true)continue;
    tout->Fill();
  }
  cout <<endl;
  tout->Write();
  fOut->Close();
  cout <<"MAPMTs data wrote\n";


}
*/

//void getSIMULATION(THeader *run);
//void getSIPM(THeader *run);
