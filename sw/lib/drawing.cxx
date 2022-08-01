#define MAXDATA 10000
#include <iostream>
#include <vector>
//#include <filesystem>
//namespace fs = std::filesystem;

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TSystem.h>

#include "definition.h"
#include "utility.h"

static TH1D *hTime;
static TH1D *hCoinc;
static TH2D *hMap;
static TH2D *hnMap;
static TH1D *hRadius;
static TH1D *hnRadius;
static vector<TH1D*> vspRadius;
static vector<TH1D*> vspTime;
static vector<TH1D*> vspPhoton;
static vector<TH1D*> vspnRadius;
static vector<TH1D*> vspnTime;
static vector<TH1D*> vspnPhoton;


void fillHisto(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pol[MAXDATA], pmt[MAXDATA], spPhoton[10], spnPhoton[10];
  double y[MAXDATA], x[MAXDATA], r[MAXDATA], nx[MAXDATA], ny[MAXDATA], nr[MAXDATA], nt[MAXDATA], spRadius[10], spTime[10], spnRadius[10], spnTime[10];
  bool coincPhoton[MAXDATA],outerPhoton[MAXDATA], goodSP[10], goodSPN[10];
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("nx",&nx);
  t->SetBranchAddress("ny",&ny);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  t->SetBranchAddress("outerPhoton",&outerPhoton);
  
  t->SetBranchAddress("spRadius",&spRadius);
  t->SetBranchAddress("spTime",&spTime);
  t->SetBranchAddress("spPhoton",&spPhoton);
  t->SetBranchAddress("goodSP",&goodSP);

  t->SetBranchAddress("spnRadius",&spnRadius);
  t->SetBranchAddress("spnTime",&spnTime);
  t->SetBranchAddress("spnPhoton",&spnPhoton);
  t->SetBranchAddress("goodSPN",&goodSPN);


  cout <<"Filling histograms\n";
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
    for(int j = 0; j < nedge; j++){
      if(pol[j]==0) hTime->Fill(nt[j]);
      if(coincPhoton[j]==true){
        hCoinc->Fill(nt[j]);
        hMap->Fill(x[j],y[j]);
        hRadius->Fill(r[j]);
        hnMap->Fill(nx[j],ny[j]);
        hnRadius->Fill(nr[j]);
      }
    }
    for(int j = 0; j < 10; j++){
      if(goodSP[j]==true){
        vspRadius[j]->Fill(spRadius[j]);
        vspTime[j]->Fill(spTime[j]);
        vspPhoton[j]->Fill(spPhoton[j]);
      }
      if(goodSPN[j]==true){
        vspnRadius[j]->Fill(spnRadius[j]);
        vspnTime[j]->Fill(spnTime[j]);
        vspnPhoton[j]->Fill(spnPhoton[j]);
        }
    }
  }
  printEnd();
  fIn->Close();
}

void displayBase(THeader *run){
  //Time distibution and coincidence peak zoom, rings before and after correction
  string out_pdf0 = Form("%soutput/plot/run%04d/displayTime.pdf[",run->suite.c_str(),run->runNum);
  string out_pdf = Form("%soutput/plot/run%04d/displayTime.pdf",run->suite.c_str(),run->runNum);
  string out_pdf1 = Form("%soutput/plot/run%04d/displayTime.pdf]",run->suite.c_str(),run->runNum);
  string out_root = Form("%soutput/plot/run%04d/displayTime.root",run->suite.c_str(),run->runNum);

  TList *save = new TList();
  save->Add(hTime);
  save->Add(hCoinc);
  save->Add(hMap);
  save->Add(hnMap);
  save->Add(hRadius);
  save->Add(hnRadius);
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  c1->Draw();
  c1->Print(out_pdf0.c_str());

  hTime->Draw();
  TLine *l1 = new TLine(run->timeMin,0,run->timeMin,hTime->GetBinContent(hTime->GetMaximumBin()));
  l1->SetLineColor(3);
  l1->Draw("same");
  TLine *l2 = new TLine(run->timeMax,0,run->timeMax,hTime->GetBinContent(hTime->GetMaximumBin()));
  l2->SetLineColor(3);
  l2->Draw("same");
  c1->Update();
  c1->Print(out_pdf.c_str());

  hCoinc->Draw();
  l1->SetY2(hCoinc->GetBinContent(hCoinc->GetMaximumBin()));
  l2->SetY2(hCoinc->GetBinContent(hCoinc->GetMaximumBin()));
  l1->Draw("same");
  l2->Draw("same");
  c1->Update();
  
  c1->Print(out_pdf.c_str());
  hMap->Draw("colz");
  c1->Update();
  c1->Print(out_pdf.c_str());
  hnMap->Draw("colz");
  c1->Print(out_pdf.c_str());
  hRadius->Draw();
  c1->Print(out_pdf.c_str());
  hnRadius->Draw();
  c1->Print(out_pdf.c_str());
  
  c1->Print(out_pdf1.c_str());

  TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
  save->Write();
  fOut->Close();
}

void displaySP(THeader *run){
  string out_pdf0 = Form("%soutput/plot/run%04d/displaySP.pdf[",run->suite.c_str(),run->runNum);
  string out_pdf = Form("%soutput/plot/run%04d/displaySP.pdf",run->suite.c_str(),run->runNum);
  string out_pdf1 = Form("%soutput/plot/run%04d/displaySP.pdf]",run->suite.c_str(),run->runNum);
  string out_root = Form("%soutput/plot/run%04d/displaySP.root",run->suite.c_str(),run->runNum);
  TList *save = new TList();
  save->Add(vspRadius[4]);
  save->Add(vspTime[4]);
  save->Add(vspPhoton[4]);
  save->Add(vspRadius[9]);
  save->Add(vspTime[9]);
  save->Add(vspPhoton[9]);
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  c1->Draw();
  c1->Print(out_pdf0.c_str());
  c1->Divide(3);
  c1->cd(1);
  vspRadius[4]->SetTitle("Single particle radius - inner ring");
  vspRadius[4]->Draw();
  c1->cd(2);
  vspTime[4]->SetTitle("Mean time of the event - inner ring");
  vspTime[4]->Draw();
  c1->cd(3);
  vspPhoton[4]->SetTitle("# photons for particle - inner ring");
  vspPhoton[4]->Draw();


  c1->Print(out_pdf.c_str());
  c1->cd(1); 
  vspRadius[9]->SetTitle("Single particle radius - outer ring");
  vspRadius[9]->Draw();
  c1->cd(2);
  vspTime[9]->SetTitle("Mean time of the event - outer ring");
  vspTime[9]->Draw();
  c1->cd(3);
  vspPhoton[9]->SetTitle("# photons for particle - outer ring");
  vspPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vspPhoton[9]->Draw();
  c1->Print(out_pdf.c_str());

  c1->Print(out_pdf1.c_str());

  TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
  save->Write();
  fOut->Close();
}

void displaySPN(THeader *run){
  string out_pdf0 = Form("%soutput/plot/run%04d/displaySPN.pdf[",run->suite.c_str(),run->runNum);
  string out_pdf = Form("%soutput/plot/run%04d/displaySPN.pdf",run->suite.c_str(),run->runNum);
  string out_pdf1 = Form("%soutput/plot/run%04d/displaySPN.pdf]",run->suite.c_str(),run->runNum);
  string out_root = Form("%soutput/plot/run%04d/displaySPN.root",run->suite.c_str(),run->runNum);
  TList *save = new TList();
  save->Add(vspnRadius[4]);
  save->Add(vspnTime[4]);
  save->Add(vspnRadius[9]);
  save->Add(vspnTime[9]);

  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  c1->Draw();
  c1->Print(out_pdf0.c_str());

  c1->Divide(3);
  c1->cd(1);
  vspnRadius[4]->SetTitle("Single particle radius - inner ring - corrected");
  vspnRadius[4]->Draw();
  c1->cd(2);
  vspnTime[4]->SetTitle("Mean time of the event - inner ring - corrected");
  vspnTime[4]->Draw();
  c1->cd(3);
  vspnPhoton[4]->SetTitle("# photons for particle - inner ring - corrected");
  vspnPhoton[4]->Draw();
  c1->Print(out_pdf.c_str());


  c1->cd(1); 
  vspnRadius[9]->SetTitle("Single particle radius - outer ring - corrected");
  vspnRadius[9]->Draw();
  c1->cd(2);
  vspnTime[9]->SetTitle("Mean time of the event - outer ring - corrected");
  vspnTime[9]->Draw();
  c1->cd(3);
  vspnPhoton[9]->SetTitle("# photons for particle - outer ring - corrected");
  vspnPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vspnPhoton[9]->Draw();
  c1->Print(out_pdf.c_str());
  c1->Print(out_pdf1.c_str());


  TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
  save->Write();
  fOut->Close();
}


void inizializePlot(THeader *run){
  //Create output directory
  string out_dir = Form("%soutput/plot/run%04d",run->suite.c_str(),run->runNum);
  if(gSystem->AccessPathName(out_dir.c_str())){
    if(std::system(Form("mkdir -p %s",out_dir.c_str()))){
      cout <<"[ERROR] Output directory not created\n";
      exit(EXIT_FAILURE);
    }
  }

  hTime = new TH1D("hTime","Time distribution for all hits",1000,0,1000);
  hCoinc = new TH1D("hCoinc","Coincidence peak",5*(int)(run->timeMax-run->timeMin),run->timeMin,run->timeMax);

  if(run->sensor=="MPPC") hMap = new TH2D("hMap","Hit position MPPC;x [mm];y [mm]",sizeof(xBinMPPC)/sizeof(*xBinMPPC)-1,xBinMPPC,sizeof(yBinMPPC)/sizeof(*yBinMPPC)-1,yBinMPPC);
  else if(run->sensor=="MAPMT") hMap = new TH2D("hMap","Hit position MAPMT;x [mm];y [mm]",sizeof(xBinMAPMT)/sizeof(*xBinMAPMT)-1,xBinMAPMT,sizeof(yBinMAPMT)/sizeof(*yBinMAPMT)-1,yBinMAPMT);
  else hMap = new TH2D("hMap","Hit position Other;x [mm];y [mm]",200,-100,100,200,-100,100);

  hnMap = new TH2D("hnMap","Corrected positions of hit;x [mm];y [mm]",200,-100,100,200,-100,100);

  hRadius = new TH1D("hRadius","Single photon radius - before corrections;r [mm]",100,0,100);
  hnRadius = new TH1D("hnRadius","Single photon radius - after corrections;r [mm]",200,0,100);

  for(int i = 0; i < 10; i++){
    TH1D *hspRadius = new TH1D(Form("hspRadius_%d",i),Form("Single particle radius - %d - before corrections;radius [mm]",i),200,0,100);
    vspRadius.push_back(hspRadius);
    TH1D *hspTime = new TH1D(Form("hspTime_%d",i),Form("Single particle radius - %d - before corrections;time [ns]",i),5*(int)(run->timeMax-run->timeMin),run->timeMin,run->timeMax);
    vspTime.push_back(hspTime);
    TH1D *hspPhoton = new TH1D(Form("hspPhoton_%d",i),Form("Single particle radius - %d - before corrections;photon [#]",i),50,0,50);
    vspPhoton.push_back(hspPhoton);
    TH1D *hspnRadius = new TH1D(Form("hspnRadius_%d",i),Form("Single particle radius - %d - after corrections;radius [mm]",i),200,0,100);
    vspnRadius.push_back(hspnRadius);
    TH1D *hspnTime = new TH1D(Form("hspnTime_%d",i),Form("Single particle radius - %d - after corrections;time [ns]",i),5*(int)(run->timeMax-run->timeMin),run->timeMin,run->timeMax);
    vspnTime.push_back(hspnTime);
    TH1D *hspnPhoton = new TH1D(Form("hspnPhoton_%d",i),Form("Single particle radius - %d - after corrections;photon [#]",i),50,0,50);
    vspnPhoton.push_back(hspnPhoton);
    cout <<i <<endl;
  }
}
