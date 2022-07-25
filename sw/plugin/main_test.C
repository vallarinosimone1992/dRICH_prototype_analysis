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
#include "../lib/readData.h"


using namespace std;

map<int,int> m1;
map<int,int>::iterator it_m1;

map<int,double> m2;
map<int,double> m3;
map<int,double>::iterator it_m2;
map<int,double>::iterator it_m3;

map<string,int> m4;
map<string,int> m5;
map<string,int>::iterator it_m4;
map<string,int>::iterator it_m5;

const double xBin[] = {-81,-77.75, -74.71875, -71.6875, -68.65625, -65.625, -62.59375, -59.5625, -56.53125, -53.5, -50.46875, -47.4375, -44.40625, -41.375, -38.34375, -35.3125, -32.28125, -29.25, -24.25, -21.21875, -18.1875, -15.15625, -12.125, -9.09375, -6.0625, -3.03125, 0, 3.03125, 6.0625, 9.09375, 12.125, 15.15625, 18.1875, 21.21875, 24.25, 29.25, 32.28125, 35.3125, 38.34375, 41.375, 44.40625, 47.4375, 50.46875, 53.5, 56.53125, 59.5625, 62.59375, 65.625, 68.65625, 71.6875, 74.71875, 77.75,81};
const double yBin[] = {-81,-77.75, -74.71875, -71.6875, -68.65625, -65.625, -62.59375, -59.5625, -56.53125, -53.5, -50.46875, -47.4375, -44.40625, -41.375, -38.34375, -35.3125, -32.28125, -29.25, -24.25, -21.21875, -18.1875, -15.15625, -12.125, -9.09375, -6.0625, -3.03125, 0, 3.03125, 6.0625, 9.09375, 12.125, 15.15625, 18.1875, 21.21875, 24.25, 29.25, 32.28125, 35.3125, 38.34375, 41.375, 44.40625, 47.4375, 50.46875, 53.5, 56.53125, 59.5625, 62.59375, 65.625, 68.65625, 71.6875, 74.71875, 77.75,81};


int main(int argc, char *argv[]){
  TApplication theApp("App",&argc,argv);
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
  //theApp.Run();
  exit(EXIT_SUCCESS);
}
