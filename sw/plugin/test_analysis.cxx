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
  gROOT->SetBatch(kTRUE);
  gStyle->SetPalette(55);
  THeader header;

  // Uncomment this to use it for quick tests
  if(argc > 1)  readHeaders(35,&header);
  else readHeaders(214,&header);
  readHeaderShort(&header);
  header.outputDir=Form("run%04d",header.runNum);


  inizializePlot(&header);
  fillHisto(&header);
  displayBase(&header);
  displaySP(&header);
  displaySPN(&header);
  displayCUT(&header);
  displayRSD(&header);

  exit(EXIT_SUCCESS);
}
