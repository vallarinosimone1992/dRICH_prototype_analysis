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
#include "utility.h"
#include "computing.h"

using namespace std;

bool correctionFit=false;
bool correctionMax=true;

double timeWalkCorrection(double t0, double t1){
    double dur = t1-t0;
    return t0+25./50.*(dur-50);
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

  double nx[MAXDATA], ny[MAXDATA], nr[MAXDATA], nrRad[MAXDATA];
  auto tNx = t->Branch("nx",&nx,"nx[nedge]/D");
  auto tNy = t->Branch("ny",&ny,"ny[nedge]/D");
  auto tNr = t->Branch("nr",&nr,"nr[nedge]/D");
  //auto tNrRad = t->Branch("nrRad",&nrRad,"nrRad[nedge]/D");

  double xNCin;
  double yNCin;
  double xNCout;
  double yNCout;
  auto txNCin = t->Branch("xNCin",&xNCin,"xNCin/D");
  auto tyNCin = t->Branch("yNCin",&yNCin,"yNCin/D");
  auto txNCout = t->Branch("xNCout",&xNCout,"xNCout/D");
  auto tyNCout = t->Branch("yNCout",&yNCout,"yNCout/D");

  cout <<"Apply position correction\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    xNCin = gxtheta*run->secondPath+run->innerCorrectionX;
    yNCin = gytheta*run->secondPath+run->innerCorrectionY;
    xNCout = gxtheta*run->firstPath+run->outerCorrectionX;
    yNCout = gytheta*run->firstPath+run->outerCorrectionY;
    for(int j = 0; j < nedge; j++){
      nx[j]=0;
      ny[j]=0;
      nr[j]=0;
      //nrRad[j]=0;
      if(coincPhoton[j]==true){
        double inPath=0, zMir;
        if(outerPhoton[j]==false){
          nx[j]=x[j]-xNCin;
          ny[j]=y[j]-yNCin;
          inPath=run->secondPath;
          zMir=run->secondMirrorPosition;
        }else{
          nx[j]=x[j]-xNCout;
          ny[j]=y[j]-yNCout;
          inPath=run->firstPath;
          zMir=run->firstMirrorPosition;
        }
        nr[j]=mmTomRad(sqrt(pow(nx[j],2)+pow(ny[j],2)),inPath,zMir);
        //nr[j]=sqrt(pow(nx[j],2)+pow(ny[j],2));
        //nrRad[j]=atan(nr[j]/path);
      }
    }
    tNx->Fill();
    tNy->Fill();
    tNr->Fill();
    //tNrRad->Fill();
    txNCin->Fill();
    tyNCin->Fill();
    txNCout->Fill();
    tyNCout->Fill();
  }
  printEnd();
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}
//////////////////////////////////////////////////////


void opticalCenterX(THeader *run){
  if(run->innerCorrectionX != 0 && run->outerCorrectionX != 0) return;
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);
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

  cout <<"Computing the x corrections\n";
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
      double dd = f0->GetParameter(1)+abs(run->innerCorrectionY);
      if(correctionFit == true && correctionMax == false) run->innerCorrectionX=f1->GetParameter(1)-dd;
      else if(correctionFit == false && correctionMax == true) run->innerCorrectionX=h1->GetBinCenter(h1->GetMaximumBin())-dd;
      else{
        cout <<"[ERROR] Undefined kind of correction on x axis. See lib/correction.cxx\n";
        exit(EXIT_FAILURE);
      }
    }
    if(run->outerCorrectionX == 0){
      TF1 *f0 = new TF1("f0","gaus(0)",-100,100);
      f0->SetParameters(100,65,1);
      h5->Fit("f0","Q","",-90,90);
      TF1 *f1 = new TF1("f1","gaus(0)",-100,100);
      f1->SetParameters(100,65,1);
      h6->Fit("f1","Q","",-90,90);
      double dd = f0->GetParameter(1)+abs(run->outerCorrectionY);
      if(correctionFit == true && correctionMax == false) run->outerCorrectionX=f1->GetParameter(1)-dd;
      else if(correctionFit == false && correctionMax == true) run->outerCorrectionX=h6->GetBinCenter(h6->GetMaximumBin())-dd;
      else{
        cout <<"[ERROR] Undefined kind of correction on x axis. See lib/correction.cxx\n";
        exit(EXIT_FAILURE);
      }
    }  
  }else{
    cout <<"[ERROR] Wrong kind of sensor from logbook\n";
    cout <<"It is " <<run->sensor.c_str() <<endl;
    exit(EXIT_FAILURE);
  }
  if(abs(run->innerCorrectionX) > 10){
    cout <<"[WARNING] Inner correction X larger than 10. Fixed to 0. Check it\n";
    run->innerCorrectionX=0;
  }
  if(abs(run->outerCorrectionX) > 10){
    cout <<"[WARNING] Outer correction X larger than 10. Fixed to 0. Check it\n";
    run->outerCorrectionX=0;
  }
  fIn->Close();
}

void opticalCenterY(THeader *run){
  if(run->innerCorrectionY != 0 && run->outerCorrectionY != 0) return;
  int runDRICH=run->runNum;
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",&run->suite[0],run->runNum);

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

  t->Draw("(spRadius[0]-spRadius[2])/2>>hIn","spPhoton[0]>0 && spPhoton[2]>0 && gxtheta<0.001 && gytheta < 0.001","goff");
  t->Draw("(spRadius[5]-spRadius[7])/2>>hOut","spPhoton[5]>0 && spPhoton[7]>0 && gxtheta<0.001 && gytheta < 0.001","goff");

  cout <<"Computing the y corrections\n";
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
  if(abs(run->innerCorrectionY) > 10){
    cout <<"[WARNING] Inner correction y larger than 10. Fixed to 0. Check it\n";
    run->innerCorrectionY=0;
  }
  if(abs(run->outerCorrectionY) > 10){
    cout <<"[WARNING] Outer correction y larger than 10. Fixed to 0. Check it\n";
    run->outerCorrectionY=0;
  }
  fIn->Close();
}

