#define MAXDATA 10000
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>

#include "definition.h"

using namespace std;

bool correctionFit=true;
bool correctionMax=false;

void newSingleParticle(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);
  TFile *fIn = new TFile(fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");
  int nedge, pmt[MAXDATA];
  double nr[MAXDATA], time[MAXDATA];
  bool coincPhoton[MAXDATA],outerPhoton[MAXDATA];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);
  
  double spnRadius[10], spnTime[10];
  int spnPhoton[10];

  auto tSpnRadius=t->Branch("spnRadius",&spnRadius,"spnRadius[10]/D");
  auto tSpnPhoton=t->Branch("spnPhoton",&spnPhoton,"spnPhoton[10]/I");
  auto tSpnTime=t->Branch("spnTime",&spnTime,"spnTime[10]/D");
 
 for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    for(int j = 0; j < 10; j++){
      spnRadius[j]=0;
      spnPhoton[j]=0;
      spnTime[j]=0;
    }
    for(int j = 0; j < nedge; j++){
      if(coincPhoton[j]==true){
        int k=0;
        if(outerPhoton[j]==true) k = 1;
        int refPMT = pmt[j]+5*k;
        int refTOT = 4+5*k;
        spnRadius[refPMT]+=nr[j];
        spnPhoton[refPMT]+=1;
        spnTime[refPMT]+=time[j];
        spnRadius[refTOT]+=nr[j];
        spnPhoton[refTOT]+=1;
        spnTime[refTOT]+=time[j];
      }
    }
    for(int j = 0; j < 10; j++){
      if(spnPhoton[j]!=0){
        spnRadius[j]/=spnPhoton[j];
        spnTime[j]/=spnPhoton[j];
        //cout <<j<<" " <<spnRadius[j]<<" " <<spnTime[j] <<" "<<spnPhoton[j] <<endl;
      }else{
        spnRadius[j]=0;
        spnTime[j]=0;
      }
    }
    //cin.get();
    tSpnRadius->Fill();
    tSpnPhoton->Fill();
    tSpnTime->Fill();
  }
  t->Write("",TObject::kOverwrite);
  fIn->Close();
  cout <<endl;
}

void positionCorrection(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);
  TFile *fIn = new TFile(fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");
  int nedge;
  float gxtheta, gytheta;
  double x[MAXDATA], y[MAXDATA];
  bool coincPhoton[MAXDATA],outerPhoton[MAXDATA];
  t->SetBranchAddress("gxtheta",&gxtheta);
  t->SetBranchAddress("gytheta",&gytheta);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);

  double nx[MAXDATA], ny[MAXDATA], nr[MAXDATA];
  auto tNx = t->Branch("nx",&nx,"nx[nedge]/D");
  auto tNy = t->Branch("ny",&ny,"ny[nedge]/D");
  auto tNr = t->Branch("nr",&nr,"nr[nedge]/D");

  double xNCin;
  double yNCin;
  double xNCout;
  double yNCout;
  auto txNCin = t->Branch("xNCin",&xNCin,"xNCin/D");
  auto tyNCin = t->Branch("yNCin",&yNCin,"yNCin/D");
  auto txNCout = t->Branch("xNCout",&xNCout,"xNCout/D");
  auto tyNCout = t->Branch("yNCout",&yNCout,"yNCout/D");

  for(int i = 0; i < t->GetEntries(); i++){
    if(i%(t->GetEntries()/10)==0)cout <<Form("\rComputing the corrected variables: %lld%% completed   ",10*(1+i)/(t->GetEntries()/10)) <<flush;
    t->GetEntry(i);
    xNCin = gxtheta*run->secondPath+run->innerCorrectionX;
    yNCin = gytheta*run->secondPath+run->innerCorrectionY;
    xNCout = gxtheta*run->firstPath+run->outerCorrectionX;
    yNCout = gytheta*run->firstPath+run->outerCorrectionY;
    for(int j = 0; j < nedge; j++){
     nx[j]=0;
     ny[j]=0;
     nr[j]=0;
      if(coincPhoton[j]==true){
        if(outerPhoton[j]==false){
          nx[j]=x[j]-xNCin;
          ny[j]=y[j]-yNCin;
        }else{
          nx[j]=x[j]-xNCout;
          ny[j]=y[j]-yNCout;
        }
        nr[j]=sqrt(pow(nx[j],2)+pow(ny[j],2));
      }
    }
    tNx->Fill();
    tNy->Fill();
    tNr->Fill();
    txNCin->Fill();
    tyNCin->Fill();
    txNCout->Fill();
    tyNCout->Fill();
  }
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}
//////////////////////////////////////////////////////


void opticalCenterX(THeader *run){
  if(run->innerCorrectionX != 0 && run->outerCorrectionX != 0) return;
  int runDRICH=run->runNum;
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],runDRICH);
  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  float gxtheta, gytheta;
  double spRadius[10], spTime[10];
  int spPhoton[10];

  t->SetBranchAddress("gxtheta",&gxtheta);
  t->SetBranchAddress("gytheta",&gytheta);
  t->SetBranchAddress("spRadius",&spRadius);
  t->SetBranchAddress("spPhoton",&spPhoton);
  t->SetBranchAddress("spTime",&spTime);

  cout <<"Computing the x corrections\n";
  TH1D *hOut = new TH1D("hOut","hOut",100,-20,20); 
  TH1D *hIn = new TH1D("hIn","hIn",100,-20,20); 
  t->Draw("(spRadius[1]-spRadius[3])/2>>hIn","spPhoton[1]>0 && spPhoton[3]>0 && gxtheta<0.001 && gytheta < 0.001","goff");
  t->Draw("(spRadius[6]-spRadius[8])/2>>hOut","spPhoton[6]>0 && spPhoton[8]>0 && gxtheta<0.001 && gytheta < 0.001","goff");

  TH1D *h0 = new TH1D("h0","h0",400,-100,100);
  TH1D *h1 = new TH1D("h1","h1",400,-100,100);
  TH1D *h5 = new TH1D("h5","h5",400,-100,100);
  TH1D *h6 = new TH1D("h6","h6",400,-100,100);

  t->Draw("spRadius[0]>>h0","spPhoton[0] > 0","goff");
  t->Draw("spRadius[1]>>h1","spPhoton[1] > 0","goff");
  t->Draw("spRadius[5]>>h5","spPhoton[5] > 0","goff");
  t->Draw("spRadius[6]>>h6","spPhoton[6] > 0","goff");

  if(run->sensor == "MAPMT"){
    if(run->innerCorrectionX == 0){
      if(correctionFit == true && correctionMax == false){
        TF1 *f = new TF1("f","gaus(0)",-20,20);
        f->SetParameters(100,0,2);
        hIn->Fit("f","Q","",-12,12);
        run->innerCorrectionX=f->GetParameter(1);
      }else if(correctionFit == false && correctionMax == true) run->innerCorrectionX=hIn->GetBinCenter(hIn->GetMaximumBin());
      else{
        cout <<"[ERROR] Undefined kind of correction on x axis. See lib/correction.cxx\n";
        exit(EXIT_FAILURE);
      }
    }
    if(run->outerCorrectionX == 0){
      if(correctionFit == true && correctionMax == false){
        TF1 *f = new TF1("f","gaus(0)",-20,20);
        f->SetParameters(100,0,2);
        hOut->Fit("f","Q","",-12,12);
        run->outerCorrectionX=f->GetParameter(1);
      }else if(correctionFit == false && correctionMax == true) run->outerCorrectionX=hOut->GetBinCenter(hOut->GetMaximumBin());
      else{
        cout <<"[ERROR] Undefined kind of correction on x axis. See lib/correction.cxx\n";
        exit(EXIT_FAILURE);
      }
    }
  }
  else if(run->sensor=="MPPC"){
    if(run->innerCorrectionX == 0){
      TF1 *f0 = new TF1("f0","gaus(0)",-100,100);
      f0->SetParameters(100,65,1);
      h0->Fit("f0","Q","",-90,90);
      TF1 *f1 = new TF1("f1","gaus(0)",-100,100);
      f1->SetParameters(100,65,1);
      h1->Fit("f1","Q","",-90,90);
      /*TF1 *f2 = new TF1("f2","gaus(0)",-100,100);
        f2->SetParameters(100,65,1);
        vh[2]->Fit("f2","Q","",-90,90);
        double dd = (f0->GetParameter(1)+abs(run->innerCorrectionY)+f1->GetParameter(1)-abs(run->innerCorrectionY))/2;*/
      double dd = f0->GetParameter(1)+abs(run->innerCorrectionY);
      if(correctionFit == true && correctionMax == false) run->innerCorrectionX=f1->GetParameter(1)-dd;
      else if(correctionFit == false && correctionMax == true) run->innerCorrectionX=h1->GetBinCenter(h1->GetMaximumBin())-dd;
      else{
        cout <<"[ERROR] Undefined kind of correction on x axis. See lib/correction.cxx\n";
        exit(EXIT_FAILURE);
      }
      cout <<"Inner "<<dd <<" " <<f1->GetParameter(1) <<" " <<h1->GetBinCenter(h1->GetMaximumBin()) <<endl;
    }
    if(run->outerCorrectionX == 0){
      TF1 *f0 = new TF1("f0","gaus(0)",-100,100);
      f0->SetParameters(100,65,1);
      h5->Fit("f0","Q","",-90,90);
      TF1 *f1 = new TF1("f1","gaus(0)",-100,100);
      f1->SetParameters(100,65,1);
      h6->Fit("f1","Q","",-90,90);
      /*TF1 *f2 = new TF1("f2","gaus(0)",-100,100);
        f2->SetParameters(100,65,1);
        vh[7]->Fit("f2","Q","",-90,90);
        double dd = (f0->GetParameter(1)+abs(run->outerCorrectionY)+f1->GetParameter(1)-abs(run->outerCorrectionY))/2;*/
      double dd = f0->GetParameter(1)+abs(run->outerCorrectionY);
      if(correctionFit == true && correctionMax == false) run->outerCorrectionX=f1->GetParameter(1)-dd;
      else if(correctionFit == false && correctionMax == true) run->outerCorrectionX=h6->GetBinCenter(h6->GetMaximumBin())-dd;
      else{
        cout <<"[ERROR] Undefined kind of correction on x axis. See lib/correction.cxx\n";
        exit(EXIT_FAILURE);
      }
      cout <<"Outer "<<dd <<" " <<f1->GetParameter(1) <<" " <<h6->GetBinCenter(h6->GetMaximumBin()) <<endl;
    }  
  }else{
    cout <<"[ERROR] Wrong kind of sensor from logbook\n";
    cout <<"It is " <<run->sensor.c_str() <<endl;
    exit(EXIT_FAILURE);
  }
  fIn->Close();
}

void opticalCenterY(THeader *run){
  if(run->innerCorrectionY != 0 && run->outerCorrectionY != 0) return;
  int runDRICH=run->runNum;
  //Find DRICH_SUITE environment variable
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],runDRICH);

  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  float gxtheta, gytheta;
  double spRadius[10], spTime[10];
  int spPhoton[10];

  t->SetBranchAddress("gxtheta",&gxtheta);
  t->SetBranchAddress("gytheta",&gytheta);
  t->SetBranchAddress("spRadius",&spRadius);
  t->SetBranchAddress("spPhoton",&spPhoton);
  t->SetBranchAddress("spTime",&spTime);

  TH1D *hOut = new TH1D("hOut","hOut",80,-20,20); 
  TH1D *hIn = new TH1D("hIn","hIn",80,-20,20); 

  cout <<"Computing the y corrections\n";
  t->Draw("(spRadius[0]-spRadius[2])/2>>hIn","spPhoton[0]>0 && spPhoton[2]>0 && gxtheta<0.001 && gytheta < 0.001","goff");
  t->Draw("(spRadius[5]-spRadius[7])/2>>hOut","spPhoton[5]>0 && spPhoton[7]>0 && gxtheta<0.001 && gytheta < 0.001","goff");

  if(run->innerCorrectionY == 0){
    if(correctionFit == true && correctionMax == false){
      TF1 *f = new TF1("f","gaus(0)",-20,20);
      f->SetParameters(100,0,2);
      hIn->Fit("f","Q","",-12,12);
      run->innerCorrectionY=f->GetParameter(1);
    }else if(correctionFit == false && correctionMax == true) run->innerCorrectionY=hIn->GetBinCenter(hIn->GetMaximumBin());
    else{
      cout <<"[ERROR] Undefined kind of correction on y axis. See lib/correction.cxx\n";
      exit(EXIT_FAILURE);
    }
  }
  if(run->outerCorrectionY == 0){
    if(correctionFit == true && correctionMax == false){ 
      TF1 *f = new TF1("f","gaus(0)",-20,20);
      f->SetParameters(100,0,2);
      hOut->Fit("f","Q","",-12,12);
      run->outerCorrectionY=f->GetParameter(1);
    }else if(correctionFit == false && correctionMax == true) run->outerCorrectionY=hOut->GetBinCenter(hOut->GetMaximumBin());
    else{
      cout <<"[ERROR] Undefined kind of correction on y axis. See lib/correction.cxx\n";
      exit(EXIT_FAILURE);
    }
  }
  fIn->Close();
}


void singleParticle(THeader *run){
  int runDRICH=run->runNum;
  //Find DRICH_SUITE environment variable
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],runDRICH);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pol[MAXDATA], time[MAXDATA], pmt[MAXDATA];
  double x[MAXDATA],y[MAXDATA],r[MAXDATA];
  bool coincPhoton[MAXDATA], outerPhoton[MAXDATA];

  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);

  double spRadius[10], spTime[10];
  int spPhoton[10];

  auto tSPRadius=t->Branch("spRadius",&spRadius,"spRadius[10]/D");
  auto tSPPhoton=t->Branch("spPhoton",&spPhoton,"spPhoton[10]/I");
  auto tSPTime=t->Branch("spTime",&spTime,"spTime[10]/D");

  cout <<"Selecting the photon\n";
  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    for(int j = 0; j < 10; j++){
      spRadius[j]=0;
      spPhoton[j]=0;
      spTime[j]=0;
    }
    for(int j = 0; j < nedge; j++){
      if(coincPhoton[j]==true){
        int k=0;
        if(outerPhoton[j]==true) k = 1;
        int refPMT = pmt[j]+5*k;
        int refTOT = 4+5*k;
        spRadius[refPMT]+=r[j];
        spPhoton[refPMT]+=1;
        spTime[refPMT]+=time[j];
        spRadius[refTOT]+=r[j];
        spPhoton[refTOT]+=1;
        spTime[refTOT]+=time[j];
      }
    }
    for(int j = 0; j < 10; j++){
      if(spPhoton[j]!=0){
        spRadius[j]/=spPhoton[j];
        spTime[j]/=spPhoton[j];
        //cout <<j<<" " <<spRadius[j]<<" " <<spTime[j] <<" "<<spPhoton[j] <<endl;
      }else{
        spRadius[j]=0;
        spTime[j]=0;
      }
    }
    //cin.get();
    tSPRadius->Fill();
    tSPPhoton->Fill();
    tSPTime->Fill();
  }
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}

