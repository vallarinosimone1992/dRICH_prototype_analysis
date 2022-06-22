#include<iostream>
#include<stdio.h>
#include<string>
#include<cstdlib>
#include<TFile.h>
#include<TString.h>
#include<TSystem.h>

#include "handleFile.h"


using namespace std;

/*void getFiles(int runDRICH, int runGEM, TTree *t){
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
    cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
    exit(EXIT_FAILURE);
  }
  TString fNamedRICH=Form("%s/DATA/dRICH_DATA/run_%04d.root",&env_var[0],runDRICH);
  TString fNameGEM=Form("%s/DATA/GEM_DATA/run_%04d_gem.root",&env_var[0],runGEM);
  cout <<Form("Analyze dRICH run: %s\n",&fNamedRICH[0]);
  if(gSystem->AccessPathName(fNamedRICH)){
    cout <<"[ERROR] dRICH run not found\n";
    return 0;
  }
  cout <<Form("Analyze GEM run: %s\n",&fNameGEM[0]);
  if(gSystem->AccessPathName(fNameGEM)){
    cout <<"[ERROR] GEM run not found\n";
    return 0;
  }
  TFile *fdRICH = new TFile(Form("%s/DATA/dRICH_DATA/run_%04d.root",&dRICHDir[0],dRICHRun),"READ");
  if(fdRICH->IsZombie()){
    cout <<"[ERROR] fdRICH is zombie! Something was wrong\n";
    return 0;
  }
  TFile *fGEM = new TFile(Form("%s/DATA/GEM_DATA/run_%04d_gem.root",&dRICHDir[0],GEMRun));
  if(fGEM->IsZombie()){
    cout <<"[ERROR] fGEM is zombie! Something was wrong\n";
    return 0;
  }



}


void getFiles(int runDRICH, int runGEM, TString *fNamedRICH, TString *fNameGEM){
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
  cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
  exit(EXIT_FAILURE);
  }
  fNamedRICH->Form("%s/DATA/dRICH_DATA/run_%04d.root",&env_var[0],runDRICH);
  fNameGEM->Form("%s/DATA/GEM_DATA/run_%04d_gem.root",&env_var[0],runGEM);
  cout <<&fNamedRICH[0] <<endl;
  cout <<&fNameGEM[0] <<endl;
  cin.get();

  }

  void getFiles(TFile *fInDRICH, int runDRICH, TFile *fInGEM, int runGEM){
  char *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
  cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
  exit(EXIT_FAILURE);
  }
  string fNameDRICH = Form("%s/DATA/dRICH_DATA/run_%04d.root",&env_var[0],runDRICH);
  cout <<"Open file " <<fNameDRICH.c_str() <<endl;
  fInDRICH = new TFile(fNameDRICH.c_str(),"READ");
  if(fInDRICH->IsZombie()){ cout <<"fInDRICH not found\n";}
  else{ cout <<"fInDRICH opened\n";}
  string fNameGEM = Form("%s/DATA/GEM_DATA/run_%04d_gem.root",&env_var[0],runGEM);
  fInGEM->Open(fNameGEM.c_str());
  auto htmp = (TH1F *) fInDRICH->Get("htd");
  htmp->Draw();
  cin.get();
  }*/
