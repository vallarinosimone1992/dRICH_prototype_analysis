#include<iostream>

#include<TSystem.h>
#include<TFile.h>

#include "utility.h"
#include "definition.h"

using namespace std;

void printProgress(double progress){
    int bar = 70;
    cout<<"[";
    int pos = bar*progress;
    for(int i = 0; i < bar; i++){
      if(i < pos) cout <<"=";
      else if(i == pos)cout  <<">";
      else cout <<" ";
    }
    cout <<"]  " <<int(progress*100) <<"%\r";
    cout.flush();
}

void printEnd(){
  int bar = 70;
  cout<<"[";
  int pos = bar;
  for(int i = 0; i < bar; i++){
    if(i < pos) cout <<"=";
    else if(i == pos)cout  <<">";
    else cout <<" ";
  }
  cout <<"] " <<int(100) <<"\n";
}

void printUsage(){
  cout <<"[WARNING] You have to run the reconstruction simply writing\n\n./reco run_number\n\nExample\n\n./reco 214\n";
  exit(EXIT_FAILURE);
}


void checkFileExistance(THeader *run){
    if(gSystem->AccessPathName(Form("%s/DATA/dRICH_DATA/run_%04d.root",&run->suite[0],run->runNum))){
    cout <<"[ERROR] dRICH run not found\n";
    exit(EXIT_FAILURE);
  }

  TFile *fdRICH = new TFile(Form("%s/DATA/dRICH_DATA/run_%04d.root",&run->suite[0],run->runNum),"READ");
  if(fdRICH->IsZombie()){
    cout <<"[ERROR] fdRICH is zombie! Something was wrong\n";
    exit(EXIT_FAILURE);
  }
  fdRICH->Close();
  if(run->runNumGEM != 0){
    if(gSystem->AccessPathName(Form("%s/DATA/GEM_DATA/run_%04d_gem.root",&run->suite[0],run->runNumGEM))){
      cout <<"[ERROR] GEM run not found\n";
      exit(EXIT_FAILURE);
    }
    TFile *fGEM = new TFile(Form("%s/DATA/GEM_DATA/run_%04d_gem.root",&run->suite[0],run->runNumGEM),"READ");
    if(fGEM->IsZombie()){
      cout <<"[ERROR] fGEM is zombie! Something was wrong\n";
      exit(EXIT_FAILURE);
    }
    fGEM->Close();
  }
}
