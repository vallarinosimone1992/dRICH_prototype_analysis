#define MAXDATA 10000
#include <iostream>
#include <vector>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TError.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TSystem.h>

#include "definition.h"
#include "utility.h"
#include "computing.h"

static TH2D *hMap;
static TH1D *hRadiusIn;
static TH1D *hnRadiusIn;
static TH1D *hTimeIn;
static TH1D *hRadiusOut;
static TH1D *hnRadiusOut;
static TH1D *hTimeOut;

static TEllipse *ringCut;
static TEllipse *ringIn;
static TEllipse *ringOut;

void fillEventDisplay(THeader *run, int ev){
  gErrorIgnoreLevel=kWarning;
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int evt, nedge, pol[MAXDATA], pmt[MAXDATA], spPhoton[10], spnPhoton[10], cutPhoton[10], fiber[MAXDATA], ch[MAXDATA], time[MAXDATA], otime[MAXDATA];
  double y[MAXDATA], x[MAXDATA], r[MAXDATA], nx[MAXDATA], ny[MAXDATA], nr[MAXDATA], nt[MAXDATA], nttw[MAXDATA], dur[MAXDATA], rsdRadius[MAXDATA], rsdTime[MAXDATA], spRadius[10], spRadiusmm[10], spTime[10], spnRadius[10], spnTime[10], cutRadius[10], cutTime[10];
  float gx0, gy0, gx1, gy1, gxa, gya, gxtheta, gytheta;
  bool trigSig[MAXDATA], goodHit[MAXDATA], goodPhoton[MAXDATA], coincPhoton[MAXDATA],externalPhoton[MAXDATA], innerPhoton[MAXDATA], outerPhoton[MAXDATA], goodSP[10], goodSPN[10], goodCUT[10];
  double trigtime, x474time, x519time, x537time;
  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("x474time",&x474time);
  t->SetBranchAddress("x519time",&x519time);
  t->SetBranchAddress("x537time",&x537time);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("fiber",&fiber);
  t->SetBranchAddress("ch",&ch);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("otime",&otime);
  t->SetBranchAddress("nx",&nx);
  t->SetBranchAddress("ny",&ny);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("nttw",&nttw);
  t->SetBranchAddress("dur",&dur);
  t->SetBranchAddress("trigSig",&trigSig);
  t->SetBranchAddress("goodHit",&goodHit);
  t->SetBranchAddress("goodPhoton",&goodPhoton);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("externalPhoton",&externalPhoton);
  //t->SetBranchAddress("innerPhoton",&innerPhoton);
  //t->SetBranchAddress("outerPhoton",&outerPhoton);
  t->SetBranchAddress("rsdRadius",&rsdRadius);
  t->SetBranchAddress("rsdTime",&rsdTime);
  t->SetBranchAddress("gx0",&gx0);
  t->SetBranchAddress("gy0",&gy0);
  t->SetBranchAddress("gx1",&gx1);
  t->SetBranchAddress("gy1",&gy1);
  t->SetBranchAddress("gxa",&gxa);
  t->SetBranchAddress("gya",&gya);
  t->SetBranchAddress("gxtheta",&gxtheta);
  t->SetBranchAddress("gytheta",&gytheta);
  t->SetBranchAddress("spRadius",&spRadius);
  t->SetBranchAddress("spRadiusmm",&spRadiusmm);
  t->SetBranchAddress("spTime",&spTime);
  t->SetBranchAddress("spPhoton",&spPhoton);
  t->SetBranchAddress("goodSP",&goodSP);
  t->SetBranchAddress("spnRadius",&spnRadius);
  t->SetBranchAddress("spnTime",&spnTime);
  t->SetBranchAddress("spnPhoton",&spnPhoton);
  t->SetBranchAddress("goodSPN",&goodSPN);
  t->SetBranchAddress("cutRadius",&cutRadius);
  t->SetBranchAddress("cutTime",&cutTime);
  t->SetBranchAddress("cutPhoton",&cutPhoton);
  t->SetBranchAddress("goodCUT",&goodCUT);

  int label=-1;
  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if(evt == ev){
      cout <<Form("Event %d found\n",ev);
      label = i;
      break;
    }
  }
  if(label < 0){
    cout <<"[ERROR] Event not found\n";
    exit(EXIT_FAILURE);
  }
  t->GetEntry(label);
  for(int i = 0; i < nedge; i++){
    if(goodPhoton[i]==true){
      hMap->Fill(x[i],y[i]);
      if(externalPhoton[i]==false){
        hRadiusIn->Fill(r[i]);
        hnRadiusIn->Fill(nr[i]);
        hTimeIn->Fill(nttw[i]);
      }else{
        hRadiusOut->Fill(r[i]);
        hnRadiusOut->Fill(nr[i]);
        hTimeOut->Fill(nttw[i]);
      }
    }
    ringCut = new TEllipse(0,0,run->geoCut);
    ringIn = new TEllipse(0,0,spRadiusmm[4]);
    ringOut = new TEllipse(0,0,spRadiusmm[9]);
  }
}

void displayEvent(THeader *run, int ev){
 
  TCanvas *c1 = new TCanvas("c1","Plots",1600,900);
  c1->Divide(3,2);
  c1->Draw();
  c1->cd(1);
  hRadiusIn->Draw();
  c1->cd(2);
  hnRadiusIn->Draw();
  c1->cd(3);
  //hTimeIn->GetXaxis()->SetRangeUser(run->timeOuMin,run->timeInMax);
  hTimeIn->Draw();
  c1->cd(4);
  hRadiusOut->Draw();
  c1->cd(5);
  hnRadiusOut->Draw();
  c1->cd(6);
  //hTimeOut->GetXaxis()->SetRangeUser(run->timeOuMin,run->timeInMax);
  hTimeOut->Draw();
  c1->Update();

  TCanvas *c0 = new TCanvas("c0","Ring",900,900);
  c0->Draw();
  hMap->SetTitle(Form("Event %d",ev));
  hMap->Draw("colz");
  ringCut->SetLineColor(kGreen);
  ringCut->SetFillStyle(0);
  ringCut->Draw("same");
  ringIn->SetLineColor(kBlue);
  ringIn->SetFillStyle(0);
  ringIn->Draw("same");
  ringOut->SetLineColor(kBlue);
  ringOut->SetFillStyle(0);
  ringOut->Draw("same");
  c0->Update();
 }

void inizializeEventDisplay(THeader *run){
  if(run->sensor=="MPPC") hMap = new TH2D("hMap","Hit position MPPC;x [mm];y [mm]",sizeof(xBinMPPC)/sizeof(*xBinMPPC)-1,xBinMPPC,sizeof(yBinMPPC)/sizeof(*yBinMPPC)-1,yBinMPPC);
  else if(run->sensor=="MAPMT") hMap = new TH2D("hMap",Form("Hit position MAPMT - run %d;x [mm];y [mm]",run->runNum),sizeof(xBinMAPMT)/sizeof(*xBinMAPMT)-1,xBinMAPMT,sizeof(yBinMAPMT)/sizeof(*yBinMAPMT)-1,yBinMAPMT);
  else hMap = new TH2D("hMap","Hit position Other;x [mm];y [mm]",180,-90,90,180,-90,90);


  hRadiusIn = new TH1D("hRadiusIn","Inner radius",100,20,70);
  hnRadiusIn = new TH1D("hnRadiusIn","Inner radius corrected",100,20,70);
  hRadiusOut = new TH1D("hRadiusOut","Outer radius",180,120,200);
  hnRadiusOut = new TH1D("hnRadiusOut","Outer radius corrected",180,120,200);
  int tMin = run->timeOuMin;
  int tMax = run->timeInMax;
  int tBin = (tMax-tMin)*4;
  hTimeIn = new TH1D("hTimeIn","Inner time",tBin,tMin,tMax);
  hTimeOut = new TH1D("hTimeOut","Outner time",tBin,tMin,tMax);
}





