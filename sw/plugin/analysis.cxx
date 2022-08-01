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
#include "../lib/drawing.h"
#include "../lib/writeHeaderText.h"


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
  gStyle->SetPalette(55);
  THeader header;
  if(argc > 1)  readHeaders(35,&header);
  else readHeaders(214,&header);
  if(argc > 1)  readHeaderShort(&header);
  else  readHeaderShort(&header);

  inizializePlot(&header);
  fillHisto(&header);
  displayBase(&header);
  displaySP(&header);
  displaySPN(&header);

  exit(EXIT_SUCCESS);
  theApp.Run();
}
