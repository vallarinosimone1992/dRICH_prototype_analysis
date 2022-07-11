#include <iostream>
#include <stdio.h>
#include <string>
#include <map>
#include <iterator>
#include <vector>

#include <TObject.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TApplication.h>

#include "../lib/fillMAPS.h"
#include "../lib/getChannel.h"
#include "../lib/MAPMTposition.h"
#include "../lib/integrate.h"


using namespace std;

map<int,int> map_GEM_rNumber;
map<int,double> map_AeroMirror_Position;
map<int,double> map_GasMirror_Position;

map<int,int>::iterator it_GEM_rNumber;
map<int,double>::iterator it_AeroMirror_Position;
map<int,double>::iterator it_GasMirror_Position;

map<string,int> map_MAPMT1;
map<string,int> map_MAPMT2;

map<string,int>::iterator it_map_MAPMT1;
map<string,int>::iterator it_map_MAPMT2;

int main(int argc, char *argv[]){
  TApplication theApp("App",&argc,argv);

  int dRICHRun;
  if(argc == 2)
    dRICHRun = atoi(argv[1]);
  else{
    cout <<"[ERROR] Missed the dRICH run number\n";
    return 0;
  }
  cout <<Form("Analysis of the dRICH run: %04d\n",dRICHRun);

  getRunNumbers(&map_GEM_rNumber,&map_AeroMirror_Position,&map_GasMirror_Position);
  getMapMAPMT(&map_MAPMT1,&map_MAPMT2);
  cout <<"Maps got\n";

  char *dRICHDir = getenv("DRICH_SUITE");
  cout <<"DRICH_SUITE PATH: "<<dRICHDir <<"\n";
  
  header runHeader;
  readHeaders(dRICHRun,&runHeader);
  int GEMRun = runHeader.runNumGEM;

  if(gSystem->AccessPathName(Form("%s/DATA/dRICH_DATA/run_%04d.root",&dRICHDir[0],dRICHRun))){
    cout <<"[ERROR] dRICH run not found\n";
    return 0;
  }
  if(GEMRun == 0){
    cout <<"[ERROR] GEM data was not taken for this run\n";
    return 0;
  }
  if(gSystem->AccessPathName(Form("%s/DATA/GEM_DATA/run_%04d_gem.root",&dRICHDir[0],GEMRun))){
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



  //theApp.Run();
  return 0;
}
