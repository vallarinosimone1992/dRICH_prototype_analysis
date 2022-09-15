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
#include "../lib/writeHeaderText.h"


using namespace std;

int main(int argc, char *argv[]){
  TApplication theApp("App",&argc,argv);
  gROOT->SetBatch(kFALSE);
  gStyle->SetPalette(55);
  vector<THeader> header;
  
  string finish;
  if(argc == 2){
    cout <<Form("Analyzing run %d\n",atoi(argv[1]));
    THeader tmpHeader;
    readHeaders(atoi(argv[1]),&tmpHeader);
    tmpHeader.outputDir=Form("run%04d",atoi(argv[1]));
    readHeaderShort(&tmpHeader);
    header.push_back(tmpHeader);
  }else if(argc >= 3){
    cout <<"Here\n";
    string output_Name=Form("%s",argv[1]);
    for(int i = 2; i < argc; i++){
      THeader tmpHeader;
      readHeaders(atoi(argv[i]),&tmpHeader);
      tmpHeader.outputDir=Form("%s",output_Name.c_str());
      readHeaderShort(&tmpHeader);
      header.push_back(tmpHeader);
    }
  }else printUsageMon();



  inizializePlot(&header[0]);
  for(int i = 0; i < header.size(); i++)  fillHisto(&header[i]);
  for(int i = 0; i < header.size(); i++)  fillHistoMon(&header[i]);
  displayMonitor(&header[0]);
  displayMonitor2(&header[0]);

  theApp.Run();
  exit(EXIT_SUCCESS);
}
