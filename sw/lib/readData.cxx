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
#include <TRandom3.h>


#include "fillMAPS.h"
#include "photoDetPosition.h"
#include "definition.h"
#include "getChannel.h"
#include "utility.h"
#include "tracking.h"

//---------------------------------------------------
void getMAPMT(THeader *runHead) {
  //---------------------------------------------------	

  int pri = 0;
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
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);

    int upEdge=0;
    for(int j = 0; j < nedge; j++){
      if(pol[j]==0)upEdge++;
    }
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
    if(DATA_2021)ntrig=1;
    if(DATA_2021)tTrigTime=3000;
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
    //if(ntrig!=1 || nX474>1 || nX519>1 || nX537>1){if(i<10)printf("Rejected event %3d (%2d %2d %2d %2d ) \n",i,ntrig,nX474,nX519,nX537);printf("\n"); continue;}


    int sedge = 0;
    for(uint j = 0; j < nedge; j++){
      if(DATA_2021==true)sedge=j;
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
      if(DATA_2021==true)sedge=j;
    }
    tNedge=sedge;
    if(wrongEvent==true){printf("\n"); continue;}
    tout->Fill();

    if(i<10 && debug){
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


void getSIMULATION(THeader *runHead){
  int pri = 0;
  
  TRandom3 rnd;
  rnd.SetSeed(123094782);

  if(runHead->sensor!="SIMULATION"){
    cout <<"[ERROR] Inconsistency between the selected data type and logbook [SENSOR]\n";
    exit(EXIT_FAILURE);
  }
  //Take the sensor and time calibration data.
  getMaps();

  TString fNamedRICH=Form("%s/DATA/SIMULATION/run_%04d.root",runHead->suite.c_str(),runHead->runNum);
  TString fOutName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",runHead->suite.c_str(),runHead->runNum);

  //Checking that files exist
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
  TTree *t = (TTree*) fdRICH->Get("htree");
  t->SetTitle("dRICH and GEM data");


  int nhits;
  double trackE[MAXDATA];
  double avg_x[MAXDATA];
  double avg_y[MAXDATA];
  double avg_z[MAXDATA];
  double avg_t[MAXDATA];
  double vx[MAXDATA];
  double vy[MAXDATA];
  double vz[MAXDATA];
  double id[MAXDATA];


  t->SetBranchAddress("nhits",&nhits);
  t->SetBranchAddress("trackE",&trackE);
  t->SetBranchAddress("avg_x",&avg_x);
  t->SetBranchAddress("avg_y",&avg_y);
  t->SetBranchAddress("avg_z",&avg_z);
  t->SetBranchAddress("avg_t",&avg_t);
  t->SetBranchAddress("vx",&vx);
  t->SetBranchAddress("vy",&vy);
  t->SetBranchAddress("vz",&vz);
  t->SetBranchAddress("id",&id);

  TFile *fOut = new TFile(fOutName,"RECREATE");
  TTree *tout = new TTree("dRICH","dRICH data");

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


  float gx0=0, gy0=0, gx1=0, gy1=0, gxa=0, gya=0, gxtheta=0, gytheta=0;
  auto tGX0=tout->Branch("gx0",&gx0,"gx0/F");
  auto tGY0=tout->Branch("gy0",&gy0,"gy0/F");
  auto tGX1=tout->Branch("gx1",&gx1,"gx1/F");
  auto tGY1=tout->Branch("gy1",&gy1,"gy1/F");
  auto tGXA=tout->Branch("gxa",&gxa,"gxa/F");
  auto tGYA=tout->Branch("gya",&gya,"gya/F");
  auto tGXtheta=tout->Branch("gxtheta",&gxtheta,"gxtheta/F");
  auto tGYtheta=tout->Branch("gytheta",&gytheta,"gytheta/F");


  double dur[MAXDATA],nttw[MAXDATA];
  bool goodHit[MAXDATA], externalPhoton[MAXDATA];;
  auto tgoodHit= tout->Branch("goodHit",&goodHit,"goodHit[nedge]/O");
  auto tdur = tout->Branch("dur",&dur,"dur[nedge]/D");
  auto tnttw = tout->Branch("nttw",&nttw,"nttw[nedge]/D");
  auto tExternalPhoton= tout->Branch("externalPhoton",&externalPhoton,"externalPhoton[nedge]/O");

  double lambda[MAXDATA];
  auto tLambda = tout->Branch("photonWavelength",&lambda,"photonWavelength[nedge]/D");

  //TO DO LIST
  //-- Ordinarli, per la ricostruzione dell'hit. Non utile nel caso delle simulazioni.
  //-- Costruire opportuno identificativo della particella con x***time
  //-- Funzione che ricostruisce slot, fiber, ch. 
  //-- Definizione di polarità, perché sia compatibile con hit reco.
  //-- Computing "tradizionale" partendo da slot-fiber-ch.
  //-- Timing


  double cmpEnergy = 900.*runHead->energyGeV;
  int cmpLUND = runHead->beamLUND;
  cout <<Form("The beam LUND is %d\nThe compare energy is %lf\n",cmpLUND,cmpEnergy);

  double trackingErrorMean=0.0;
  double trackingErrorSigma=0.0;
  if(APPLY_SIMULATION_TRACKING_ERROR) trackingErrorSigma=0.4;
  
  //std::default_random_engine trackingErrorGenerator;
  //std::normal_distribution<double> trackingErrorDistribution(trackingErrorMean,trackingErrorSigma);

  cout <<"Taking the simulation data\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    if(debug)cout <<Form("");
    if(debug)cout <<Form("");
    t->GetEntry(i);
    tnedge=0; 
    bool flagUpGEM=false, flagDownGEM=false;
    for(int j = 0; j < nhits; j++){
      if((int)id[j]==cmpLUND && abs(avg_t[j] - UPSTREAM_GEM_TIME_SIM) < .2 && abs(avg_z[j] - UPSTREAM_GEM_Z_SIM) < 1 && trackE[j] > cmpEnergy){
        gx0=avg_x[j]+rnd.Gaus(trackingErrorMean,trackingErrorSigma);
        gy0=avg_y[j]+rnd.Gaus(trackingErrorMean,trackingErrorSigma);
        flagUpGEM=true;
      }
      if((int)id[j]==cmpLUND && abs(avg_t[j] - DOWNSTREAM_GEM_TIME_SIM) < .2 && abs(avg_z[j] - DOWNSTREAM_GEM_Z_SIM) < 1 && trackE[j] > cmpEnergy ){
        gx1=avg_x[j]+rnd.Gaus(trackingErrorMean,trackingErrorSigma);
        gy1=avg_y[j]+rnd.Gaus(trackingErrorMean,trackingErrorSigma);
        flagDownGEM=true;
      }
      if(flagUpGEM==true && flagDownGEM==true)break;
    }
    if(flagUpGEM==false || flagDownGEM==false)continue;
    if(debug)cout <<Form("GEM: %f %f %lf %lf %lf %lf %lf %lf %lf %lf\n",runHead->UpGEMz,runHead->DnGEMz,gx0,gy0,gx1,gy1,gxa,gya,gxtheta,gytheta);
    for(int j = 0; j < nhits; j++){
      if(abs(avg_z[j]-SIMULATION_DETECTOR_Z)<.1 && id[j]==0 && (abs(avg_t[j]-SIMULATION_GAS_PHOTON_EXPECTED_TIME)<0.5 || (abs(avg_t[j]-SIMULATION_AEROGEL_PHOTON_EXPECTED_TIME)<0.5))){
        //It's a photon!
        if(simulationPixel(avg_x[j],avg_y[j],&tpmt[tnedge],&tx[tnedge],&ty[tnedge],&trMM[tnedge])==false)continue;
        double path = SIMULATION_GAS_PATH;
        externalPhoton[tnedge]=false;
        if(trMM[tnedge] > runHead->geoCut){
          path= SIMULATION_AERO_PATH;
          externalPhoton[tnedge]=true;
        }
        tr[tnedge]=1000*atan(trMM[tnedge]/path);
        ttime[tnedge]=avg_t[j];
        totime[tnedge]=avg_t[j];
        tnt[tnedge]=avg_t[j];
        if(debug)cout <<Form("Hit: %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",tnedge,tpmt[j],avg_x[j],tx[tnedge],avg_y[j],ty[tnedge],trMM[tnedge],path,tr[tnedge],tnt[tnedge]);
        tpol[tnedge]=0;
        tslot[tnedge]=-1;
        tfiber[tnedge]=-1;
        tch[tnedge]=-1;
        tboard[tnedge]=-1;
        tchip[tnedge]=-1;
        ttrigSig[tnedge]=false;
        if(APPLY_QUANTUM_EFFICIENCY==false) goodHit[tnedge]=true; //Can introduce here the quantum efficiency
        else goodHit[tnedge]=applyQuantumEfficiency(trackE[j]);
        dur[tnedge]=1.1*CUT_MIN_DUR;
        nttw[tnedge]=tnt[tnedge];
        lambda[tnedge]=1240./trackE[j]/1e6;
        tnedge++;
      } 
    }

    AERO_computing(runHead,&gxa,&gya,&gxtheta,&gytheta,gx0,gy0,gx1,gy1);


    tout->Fill();
  }
  printEnd();
  tout->Write();
  fOut->Close();
}

