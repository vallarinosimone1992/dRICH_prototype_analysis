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

//---------------------------------------------------
void getMAPMT(THeader *runHead) {
//---------------------------------------------------	

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

    int evt, tPol[MAXDATA], tSlot[MAXDATA], tFiber[MAXDATA], tCh[MAXDATA] ,tTime[MAXDATA], board[MAXDATA], chip[MAXDATA], oTime[MAXDATA];
    uint tNedge;
    double tTrigTime, tX474Time, tX519Time, tX537Time, trigger;
    double x[MAXDATA], y[MAXDATA], radius[MAXDATA],nt[MAXDATA];
    int pmt[MAXDATA], anode[MAXDATA];
    bool trigSig[MAXDATA];

    auto evento=tout->Branch("evt",&evt,"evt/I");
    (void)evento;
    auto ttrigtime=tout->Branch("trigtime",&tTrigTime,"trigtime/D");
    (void)ttrigtime;
    auto tx474time=tout->Branch("x474time",&tX474Time,"x474time/D");
    (void)tx474time;
    auto tx519time=tout->Branch("x519time",&tX519Time,"x519time/D");
    (void)tx519time;
    auto tx537time=tout->Branch("x537time",&tX537Time,"x537time/D");
    (void)tx537time;
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
    auto totime=tout->Branch("origTime",&oTime,"original_time[nedge]/I");
    (void)totime;
    auto tpmt=tout->Branch("pmt",&pmt,"pmt[nedge]/I");
    (void)tpmt;
    auto tanode=tout->Branch("anode",&anode,"anode[nedge]/I");
    (void)tanode;
    auto tx=tout->Branch("x",&x,"x[nedge]/D");
    (void)tx;
    auto ty=tout->Branch("y",&y,"y[nedge]/D");
    (void)ty;
    auto tradius=tout->Branch("radius",&radius,"radius[nedge]/D");
    (void)tradius;
    auto tnt=tout->Branch("nt",&nt,"nt[nedge]/D");
    (void)tnt;
    auto ttrigSig=tout->Branch("trigSig",&trigSig,"trigSig[nedge]/O");
    (void)ttrigSig;

    cout <<"Taking the dRICH prototype data\n";
    for(int i = 0; i < t->GetEntries(); i++){
        if(SHOW_PROGRESS==true && i%100==0)printProgress((double)i/t->GetEntries());
        t->GetEntry(i);
       
        int upEdge=0;
        for(int j = 0; j < nedge; j++){
            if(pol[j]==0)upEdge++;
        }
        /*if(i<10){
            printf("\n \n PRIMA eventi %3d %6d\n",i,evtDRICH);
            for(int j = 0; j < nedge; j++)
	        printf("e     %4d  %3d (%3d %3d %d) %6d \n",j,pol[j],slot[j],fiber[j],chM[j],time[j] );
        }*/
        for(int j = 0; j < nedge; j++){
            for(int k = 0; k < nedge; k++){
                if(j == k)continue;
                if((time[j] < time[k] && pol[j]==pol[k]) || pol[j] < pol[k]){
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
        /*if(i<10){
            printf("DOPO \n \n");
            for(int j = 0; j < nedge; j++)
	        printf("e     %4d  %3d (%3d %3d %d) %6d \n",j,pol[j],slot[j],fiber[j],chM[j],time[j] );
        }*/

        evt=(int)evtDRICH;
        bool wrongEvent=false;
        int ntrig=0;
        int nX474=0;
        int nX519=0;
        int nX537=0;
        tTrigTime=0;
        tX474Time=0;
        tX519Time=0;
        tX537Time=0;
        for(int j = 0; j < nedge; j++){
            if(fiber[j]==11 && chM[j]==27 && pol[j]==0){
                ntrig++;
                tTrigTime=time[j];
                if(i<10)printf("Found     trigger  %3d time %6d \n",j,time[j]);
            }
            if(fiber[j]==11 && chM[j]==4 && pol[j]==0){
                nX474++;
                tX474Time=time[j];
                if(i<10)printf("Found     XCET_474 %3d time %6d \n",j,time[j]);
            }
            if(fiber[j]==11 && chM[j]==2 && pol[j]==0){
                nX519++;
                tX519Time=time[j];
                if(i<10)printf("Found     XCET_519 %3d time %6d \n",j,time[j]);
            }
            if(fiber[j]==11 && chM[j]==32 && pol[j]==0){
                nX537++;
                tX537Time=time[j];
                if(i<10)printf("Found     XCET_537 %3d time %6d \n",j,time[j]);
            }
        }
        if(ntrig!=1 || nX474>1 || nX519>1 || nX537>1){if(i<10)printf("Rejected event %3d (%2d %2d %2d %2d ) \n",i,ntrig,nX474,nX519,nX537);continue;}

	int sedge = 0;
        for(uint j = 0; j < nedge; j++){
            tPol[sedge]=pol[j];
            tSlot[sedge]=slot[j];
            tFiber[sedge]=fiber[j];
            tCh[sedge]=chM[j];
            chip[sedge]=getMarocChip(chM[j]);
            //if(fiber[j]==11 && chip[sedge]==0)continue;

            board[sedge]=getMarocBoard(fiber[j],runHead);
            upstreamMaroc(fiber[j], runHead);
            anode[sedge]= getMAPMT_anode(fiber[j],chM[j],board[sedge],chip[sedge],runHead->upstreamBoard);
            oTime[sedge]=time[j];
            if(fiber[j]!=11 || chip[sedge]!=0){
                tTime[sedge]=time[j] - tTrigTime + runHead->GlobalTimeOff;
                trigSig[sedge]=false;
            }else {
                tTime[sedge]=time[j] + runHead->GlobalTimeOff;
                trigSig[sedge]=true;
            }
            pmt[sedge]=FiberToPhDet(fiber[j],&(runHead->fiberRef)[0]);
            MAPMTposition(anode[sedge],pmt[sedge],&x[j],&y[j],&radius[j]);
            nt[sedge]=tTime[sedge];
            if(trigSig[sedge]==false)nt[sedge]=timeCalibrationMAPMT(tTime[sedge],anode[sedge],pmt[sedge]);
	    sedge++;
        }
        tNedge=sedge;
        if(wrongEvent==true)continue;
        tout->Fill();

        if(i<10){
            printf("NEW tree event %3d %3d nedge %5d trigger %7.2f Cher %7.2f %7.2f %7.2f \n",i,evt,tNedge,tTrigTime,tX474Time,tX519Time,tX537Time);
            for(int j = 0; j < tNedge; j++)
	        printf("e     %4d  %3d (%3d %3d %3d --> %2d) chip %3d anode %6d T %6d  %7.2f\n",j,tPol[j],tSlot[j],tFiber[j],tCh[j],
				      trigSig[j],chip[j],anode[j],tTime[j],nt[j]);
        }

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
