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
#include <TROOT.h>

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
#include "../lib/eventDisplay.h"
#include "../lib/writeHeaderText.h"


using namespace std;

int main(int argc, char *argv[]){
  TApplication theApp("App",&argc,argv);
  //gROOT->SetBatch(kTRUE);
  gStyle->SetPalette(55);
  vector<THeader> header;

  if(argc != 3){
    cout <<"[ERROR] the evente display require two arguments: run_number and event_number.\n";
    exit(EXIT_FAILURE);
  }else{
    cout <<Form("Analyzing event %d of run %d\n",atoi(argv[2]),atoi(argv[1]));
    THeader tmpHeader;
    readHeaders(atoi(argv[1]),&tmpHeader);
    tmpHeader.outputDir=Form("run%04d_event%06d",atoi(argv[1]),atoi(argv[2]));
    readHeaderShort(&tmpHeader);
    header.push_back(tmpHeader);
  }

  

  //Inizialize the event plot
  //inizializeEvent(&header[0]);
  //Fill the single event plot
  //fillEvent(&header[0]);
  //Show the single event plot
  inizializeEventDisplay(&header[0]);
  fillEventDisplay(&header[0], atoi(argv[2]));
  displayEvent(&header[0], atoi(argv[2]));

  theApp.Run();
  exit(EXIT_SUCCESS);
}
