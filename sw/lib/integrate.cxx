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
#include "getChannel.h"
#include "MAPMTposition.h"

using namespace std;

//void step1(TTree *t, TTree *tGEM, TTree *tout){
void TTreeIntegration(int runDRICH, int runGEM){
  const char  *tmp = getenv("DRICH_SUITE");
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

  //SOME USEFUL DECLARATION, BEFORE FINISH TO IMPLEMENT THE HEADER READING
  string phDet="MAPMT";


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
  cout <<"GEM entries setted\n";

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
        if(fiber[j] == 4 || fiber[j] == 5 || fiber[j] == 6 || fiber[j] == 7 ){ //To be defined
                                                                               //channel[j]=MAPMT_channel_two(fiber[j],chM[j]);
          channel[j]= 10;
        }else{
          //channel[j]=MAPMT_channel_three(fiber[j],chM[j]);
          channel[j]=20;
        }
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
  cout <<"Entries before save: "<<tout->GetEntries() <<endl;
  fOut->Close();
  cout <<"End of the function\n";
}

/*int MAPMT_channel_three(int fpga, int mChannel){
  string label;
  int channel=300;
  if(fpga == 4 || fpga == 5 || fpga == 6 || fpga == 7){
  if(mChannel > 63 && mChannel < 128){
  label=Form("IN1_%02d",mChannel-64);
  channel=map_MAPMT1.at(label.c_str());
  }
  if(mChannel>127){
  label=Form("IN3_%02d",mChannel-128);
  channel=map_MAPMT1.at(label.c_str());
  }
  }
  if(fpga == 8 || fpga == 9 || fpga == 10 || fpga == 11){
  if(mChannel > 63 && mChannel < 128){
  label=Form("IN1_%02d",mChannel-64);
  channel=map_MAPMT2.at(label.c_str());
  }
  if(mChannel>127){
  label=Form("IN3_%02d",mChannel-128);
  channel=map_MAPMT2.at(label.c_str());
  }
  }
  return channel;
  }

  int MAPMT_channel_two(int fpga, int mChannel){
  string label;
  int channel=300;
  if(fpga == 4 || fpga == 5 || fpga == 6 || fpga == 7){
  if(mChannel < 64){
  label=Form("IN1_%02d",mChannel);
  channel=map_MAPMT1.at(label.c_str());
  }
  if(mChannel>127){
  label=Form("IN3_%02d",mChannel-128);
  channel=map_MAPMT1.at(label.c_str());
  }
  }
  if(fpga == 8 || fpga == 9 || fpga == 10 || fpga == 11){
  if(mChannel < 64){
  label=Form("IN1_%02d",mChannel);
  channel=map_MAPMT2.at(label.c_str());
  }
  if(mChannel>127){
  label=Form("IN3_%02d",mChannel-128);
  channel=map_MAPMT2.at(label.c_str());
  }
  }
  return channel;
  }


  int MPPC_channel_two(int fpga, int mChannel){
  string label;
  int channel=300;
  if(fpga == 0 || fpga == 4 || fpga == 5){
  if(mChannel < 64){
  label=Form("IN1_%d",mChannel);
  channel=map_MPPC1.at(label.c_str());

  }
  if(mChannel>127){
  label=Form("IN3_%d",mChannel-128);
  channel=map_MPPC1.at(label.c_str());
  }
  }
  if(fpga == 1 || fpga == 6 || fpga == 7){
  if(mChannel < 64){
  label=Form("IN1_%d",mChannel);
  channel=map_MPPC2.at(label.c_str());
}
if(mChannel>127){
  label=Form("IN3_%d",mChannel-128);
  channel=map_MPPC2.at(label.c_str());
}
}
return channel;
}
*/
