#define MAXDATA 10000
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
#include <TEllipse.h>
#include <TAxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TError.h>
#include <TMath.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TRandom3.h>
#include <TDirectory.h>

#include "../lib/definition.h"
#include "../lib/fillMAPS.h"
#include "../lib/getChannel.h"
#include "../lib/photoDetPosition.h"
#include "../lib/selection.h"
#include "../lib/correction.h"
#include "../lib/integrate.h"
#include "../lib/computing.h"
#include "../lib/readData.h"


using namespace std;


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
  
  header info;
  readHeaders(dRICHRun,info);
  if(gSystem->AccessPathName(Form("%s/DATA/dRICH_DATA/run_%04d.root",info.suite.c_str(),info.runNum))){
    cout <<Form("[ERROR] dRICH run %d not found\n",info->runNum);
    return 0;
  }
  if(info.runNumGEM != 0){
    if(gSystem->AccessPathName(Form("%s/DATA/GEM_DATA/run_%04d_gem.root",info.suite.c_str(),info.runNumGEM))){
      cout <<Form("[ERROR] GEM run %d not found\n",info.runNumGEM);
      return 0;
    }
  }



  TCanvas *c = new TCanvas();
  c->Draw();
  TH2D *h = new TH2D("h","Rings; x[mm]; y[mm]",sizeof(xBin)/sizeof(*xBin)-1,xBin,sizeof(yBin)/sizeof(*yBin)-1,yBin);
  getRunNumbers(&m1, &m2, &m3);
  getMapMAPMT(&m4,&m5);

  THeader runHeader;
  if(argc > 1)  readHeaders(35,&runHeader);
  else readHeaders(214,&runHeader);
  cout <<" number from header: "<<runHeader.runNum <<endl;
  cout <<Form("Reading some info from run %d header file\n",runHeader.runNum);
  cout <<Form("The GEM run is %d\n",runHeader.runNumGEM);
  cout <<Form("There was a %s beam of %d GeV\n",(runHeader.beam).c_str(),runHeader.energyGeV);
  cout <<Form("Sensors were %s\n",(runHeader.sensor).c_str());
  cout <<Form("The paths are %lf %lf \n", runHeader.firstPath, runHeader.secondPath); 
  cout <<Form("DRICH_SUITE = %s", runHeader.suite.c_str());
  //return 0;

  if(runHeader.sensor=="MAPMT")getMAPMT(&runHeader);
  if(runHeader.sensor=="MPPC")getMPPC(&runHeader);



  //TTreeIntegration(runHeader.runNum,runHeader.runNumGEM);
  TTreeIntegration(&runHeader);
  cout <<"Integration done\n";

  findTimeCoincidence(&runHeader);
  cout <<"Time coincidence window extremes: " <<runHeader.timeMin <<" " <<runHeader.timeMax <<endl;
  selectPhotons(&runHeader);
  cout <<"Photon selected\n";
  singleParticle(&runHeader);
  cout <<"Single particle quantities computed\n";
  //IMPORTANT! The Y axis correction must be computed befor the X axis correction.
  opticalCenterY(&runHeader);
  opticalCenterX(&runHeader);
  cout <<"Y correction in and out: " <<runHeader.innerCorrectionY <<" " <<runHeader.outerCorrectionY <<endl; 
  cout <<"X correction in and out: " <<runHeader.innerCorrectionX <<" " <<runHeader.outerCorrectionX <<endl;
  positionCorrection(&runHeader);
  newSingleParticle(&runHeader);
  computeRMS(&runHeader);
  rmsCutSelection(&runHeader);
  computeCutSingleParticle(&runHeader);
  //theApp.Run();
  exit(EXIT_SUCCESS);
}
