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
#include "../lib/utility.h"
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
  //TApplication theApp("App",&argc,argv);

  int run;
  if(argc == 2)
    run = atoi(argv[1]);
  else printUsage();
  cout <<Form("Analysis of the dRICH run: %04d\n",run);
  
  THeader header;//It is defined in ../lib/definition.h
  readHeaders(run,&header); //fillMaps.cxx
  checkFileExistance(&header);
  cout <<Form("Reading some info from run %d header file\n",header.runNum);
  cout <<Form("The GEM run is %d\n",header.runNumGEM);
  cout <<Form("There was a %s beam of %d GeV\n",(header.beam).c_str(),header.energyGeV);
  cout <<Form("Sensors were %s\n",(header.sensor).c_str());
  cout <<Form("The paths are %lf %lf \n", header.firstPath, header.secondPath); 
  cout <<Form("DRICH_SUITE = %s", header.suite.c_str());
  //return 0;

  //Convert dRICH raw data in position and calibrated time for MAPMT and MPPC
  if(header.sensor=="MAPMT")getMAPMT(&header); //readData.cxx
  else if(header.sensor=="MPPC")getMPPC(&header); //readData.cxx
  else{
    cout <<"[ERROR] sensor not recognized\n";
    exit(EXIT_FAILURE);
  }

  //Integration betweed dRICH data and tracking data.
  if(header.runNumGEM!=0)TTreeIntegration(&header);//integration.cxx
  else noGEM_Integration(&header);


  //Find the time coincidence window extreme
  findTimeCoincidence(&header); //selection.cxx

  //Select the photons inside the time coincidence window and distinguish between inner and outer (waiting a better way to do this, based also on time).
  selectPhotons(&header); //selection.cxx

  //Compute the mean quantities for the single particle.
  singleParticle(&header); //computing.cxx

  //IMPORTANT! The Y axis correction must be computed befor the X axis correction.
  //Compute the correction factor base on the optical center.
  opticalCenterY(&header); //correction.cxx
  opticalCenterX(&header); //correction.cxx

  //Apply the position correction.
  positionCorrection(&header); //correction.cxx
                               //Compute the mean quantities for the single particle after the position correction
  newSingleParticle(&header); //computing.cxx
                              //Calculate the time and radius rms, and the residue for each photon.
  computeRMS(&header); //computing.cxx
                       //Apply the cut based on rms
  rmsCutSelection(&header); //selection.cxx
                            //Compute the mean quantities for the single particle after the rms cut application
  computeCutSingleParticle(&header); //computing.cxx

  writeHeader(&header);
  exit(EXIT_SUCCESS);
}
