#define MAXDATA 10000
#include <iostream>
#include <vector>
//#include <filesystem>
//namespace fs = std::filesystem;

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TError.h>
#include <TStyle.h>
#include <TSystem.h>

#include "definition.h"
#include "utility.h"
#include "computing.h"

static TH1D *hTime;
static TH1D *hCoinc;
static TH2D *hMap;
static TH2D *hnMap;
static TH1D *hRadius;
static TH1D *hnRadius;

static vector<TH1D*> vrsdRadius;
static vector<TH1D*> vrsdTime;

static vector<TH1D*> vspRadius;
static vector<TH1D*> vspTime;
static vector<TH1D*> vspPhoton;
static vector<TH1D*> vspnRadius;
static vector<TH1D*> vspnTime;
static vector<TH1D*> vspnPhoton;
static vector<TH1D*> vcutRadius;
static vector<TH1D*> vcutTime;
static vector<TH1D*> vcutPhoton;

static vector<TH1D*> vspSigPhoGas;
static vector<TH1D*> vspSigPhoAero;
static vector<TH1D*> vspnSigPhoGas;
static vector<TH1D*> vspnSigPhoAero;
static vector<TH1D*> vcutSigPhoGas;
static vector<TH1D*> vcutSigPhoAero;
static TH1D *hspSigVsPhoGas;
static TH1D *hspnSigVsPhoGas;
static TH1D *hcutSigVsPhoGas;
static TH1D *hspSigVsPhoAero;
static TH1D *hspnSigVsPhoAero;
static TH1D *hcutSigVsPhoAero;

static int minPhoGas=2;
static int maxPhoGas=50;
static int minPhoAero=2;
static int maxPhoAero=15;



void fillHisto(THeader *run){
  gErrorIgnoreLevel=kWarning;
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int nedge, pol[MAXDATA], pmt[MAXDATA], spPhoton[10], spnPhoton[10], cutPhoton[10];
  double y[MAXDATA], x[MAXDATA], r[MAXDATA], nx[MAXDATA], ny[MAXDATA], nr[MAXDATA], nt[MAXDATA], rsdRadius[MAXDATA], rsdTime[MAXDATA], spRadius[10], spTime[10], spnRadius[10], spnTime[10], cutRadius[10], cutTime[10];
  bool coincPhoton[MAXDATA],outerPhoton[MAXDATA], goodSP[10], goodSPN[10], goodCUT[10];
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
  t->SetBranchAddress("rsdRadius",&rsdRadius);
  t->SetBranchAddress("rsdTime",&rsdTime);
  
  
  t->SetBranchAddress("spRadius",&spRadius);
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


  cout <<Form("Filling histograms for run %d\n",run->runNum);;
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
  }
  printEnd();
  fIn->Close();
}

void displayThetaAna(THeader *run){
  gErrorIgnoreLevel=kWarning;
  //Time distibution and coincidence peak zoom, rings before and after correction
  string out_pdf0 = Form("%soutput/plot/%s/displayMonitor.pdf[",run->suite.c_str(),run->outputDir.c_str());
  string out_pdf = Form("%soutput/plot/%s/displayMonitor.pdf",run->suite.c_str(),run->outputDir.c_str());
  string out_pdf1 = Form("%soutput/plot/%s/displayMonitor.pdf]",run->suite.c_str(),run->outputDir.c_str());
  string out_root = Form("%soutput/plot/%s/displayMonitor.root",run->suite.c_str(),run->outputDir.c_str());

  TList *save = new TList();
  save->Add(hTime);

  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  c1->Divide(3,2);
  c1->Draw();
  c1->Print(out_pdf0.c_str());

  c1->Update();
  c1->Print(out_pdf.c_str());


  c1->Print(out_pdf1.c_str());
  c1->Close();

  TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
  save->Write();
  fOut->Close();
}


void inizializePlot(THeader *run){
  //Create output directory
  string out_dir = Form("%soutput/plot/%s",run->suite.c_str(),run->outputDir.c_str());
  if(gSystem->AccessPathName(out_dir.c_str())){
    if(std::system(Form("mkdir -p %s",out_dir.c_str()))){
      cout <<"[ERROR] Output directory not created\n";
      exit(EXIT_FAILURE);
    }
  }

  hRadius = new TH1D("hRadius","Single photon radius - before corrections;r [mRad]",400,0,200);
  hnRadius = new TH1D("hnRadius","Single photon radius - after corrections;r [mRad]",400,0,200);


  for(int i = 0; i < thetaBin; i++){
    TH1D *hcutRad_Theta= new TH1D(Form("hcutRad_Theta_bin_%02d",i),Form("hcut Radius - #{theta} bin %02d",i),400,0,200);
    vcutRad_Theta.push_back(hcutRad_Theta);
    TH1D *hcutSig_Theta= new TH1D(Form("hcutSig_Theta_bin_%02d",i),Form("hcut #sigma - #{theta} bin %02d",i),400,0,200);
    vcutSig_Theta.push_back(hcutSig_Theta);
  }
    TH1D *hcutSigPhoAero= new TH1D(Form("hcutSigPhoAero_%02d",i),Form("hcutSigPhoAero_%02d",i),400,0,200);
    vcutSigPhoAero.push_back(hcutSigPhoAero);
  }
  hspSigVsPhoGas = new TH1D("hspSigVsPhoGas","#sigma_{r} vs photon per particle - Gas;Photon/particle [#];#sigma_{r} [mRad]",maxPhoGas+1,-0.5,maxPhoGas+0.5);
  hspnSigVsPhoGas= new TH1D("hsnpSigVsPhoGas","#sigma_{r} vs photon per particle - Gas;Photon/particle [#];#sigma_{r} [mRad]",maxPhoGas+1,-0.5,maxPhoGas+0.5);
  hcutSigVsPhoGas = new TH1D("hcutSigVsPhoGas","#sigma_{r} vs photon per particle - Gas;Photon/particle [#];#sigma_{r} [mRad]",maxPhoGas+1,-0.5,maxPhoGas+0.5);
  hspSigVsPhoAero = new TH1D("hspSigVsPhoAero","#sigma_{r} vs photon per particle - Aerogel;Photon/particle [#];#sigma_{r} [mRad]",maxPhoAero+1,-0.5,maxPhoAero+0.5);
  hspnSigVsPhoAero= new TH1D("hsnpSigVsPhoAero","#sigma_{r} vs photon per particle - Aerogel;Photon/particle [#];#sigma_{r} [mRad]",maxPhoAero+1,-0.5,maxPhoAero+0.5);
  hcutSigVsPhoAero = new TH1D("hcutSigVsPhoAero","#sigma_{r} vs photon per particle - Aerogel;Photon/particle [#];#sigma_{r} [mRad]",maxPhoAero+1,-0.5,maxPhoAero+0.5);
}
