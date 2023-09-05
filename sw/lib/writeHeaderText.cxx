#define MAXLINEA 1000
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>


#include <TTree.h>
#include <TFile.h>
//#include <TSystem.h>

#include "writeHeaderText.h"

#include "utility.h"
#include "definition.h"

void writeHeaderTTree(THeader *run){
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"UPDATE");
  TTree *t = (TTree*) fIn->Get("dRICH");

  //THeader *heade;
  auto theader = t->Branch("header",&run,"runNum/I:firstMirrorPosition/F:secondMirrorPosition/F");
  for(int i = 0; i <<t->GetEntries(); t++){
    theader->Fill();
  }
  t->Write("",TObject::kOverwrite);
  fIn->Close();
}

void readHeaderShort(THeader *run){
  FILE *f;
  string fname=Form("%s/output/header/run%04d_header_data.txt",run->suite.c_str(),run->runNum);
  f=fopen(fname.c_str(),"r");
  if(f==NULL){
    cout <<"[ERROR] Header recap not found. Do you executed the writeHeaderShort in your reco?\n";
    exit(EXIT_FAILURE);
  }
  char line[MAXLINEA];
  char tmp[MAXLINEA];

  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->runNum)); 
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->day=tmp;}
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->startTime=tmp;}
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->endTime=tmp;}
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->beam=tmp;}
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->energyGeV));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->expEvents));
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->sensor=tmp;}
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->firstMirrorPosition));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->secondMirrorPosition));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->temperature));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->powerHV));
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->trigger=tmp;}
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->runType=tmp;}
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->runNumGEM));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->pedestalGEM));
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->setupFile=tmp;}
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->note=tmp;}
  if(fgets(line,MAXLINEA,f)){sscanf(line,"%s",tmp);
    run->suite=tmp;}
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->fiberRef[0]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->fiberRef[1]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->fiberRef[2]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->fiberRef[3]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->fiberRef[4]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->fiberRef[5]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->fiberRef[6]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->fiberRef[7]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->marocBoard[0]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->marocBoard[1]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->marocBoard[2]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->marocBoard[3]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->marocBoard[4]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->marocBoard[5]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->marocBoard[6]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->marocBoard[7]));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->upstreamBoard));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->firstPath));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->secondPath));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->UpGEMz));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->DnGEMz));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->zAerogel));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->geoCut));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->radCut));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->cutRadiusInRMS));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->cutTimeInRMS));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->cutRadiusOutRMS));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->cutTimeOutRMS));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->innerCorrectionX));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->innerCorrectionY));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->outerCorrectionX));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->outerCorrectionY));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->UpGEMxRunOff));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->UpGEMyRunOff));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->DnGEMxRunOff));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%f",&(run->DnGEMyRunOff));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->timeInMin));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->timeInMax));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->timeOuMin));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->timeOuMax));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->px474));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->px519));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->px537));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->beamChLogic));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->lookbackDAQ));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->aerogelRefractiveIndex));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->MaxHitLength));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->GlobalTimeOff));
  
  fclose(f);
}




void writeHeaderShort(THeader *run){
  FILE *f;
  string fname=Form("%s/output/header/run%04d_header_data.txt",run->suite.c_str(),run->runNum);
  f=fopen(fname.c_str(),"w");

  fprintf(f,"%d\n",(run->runNum)); 
  fprintf(f,"%s\n",run->day.c_str());
  fprintf(f,"%s\n",run->startTime.c_str());
  fprintf(f,"%s\n",run->endTime.c_str());
  fprintf(f,"%s\n",run->beam.c_str());
  fprintf(f,"%lf\n",(run->energyGeV));
  fprintf(f,"%d\n",(run->expEvents));
  fprintf(f,"%s\n",run->sensor.c_str());
  fprintf(f,"%f\n",(run->firstMirrorPosition));
  fprintf(f,"%f\n",(run->secondMirrorPosition));
  fprintf(f,"%f\n",(run->temperature));
  fprintf(f,"%d\n",(run->powerHV));
  fprintf(f,"%s\n",run->trigger.c_str());
  fprintf(f,"%s\n",run->runType.c_str());
  fprintf(f,"%d\n",(run->runNumGEM));
  fprintf(f,"%d\n",(run->pedestalGEM));
  fprintf(f,"%s\n",run->setupFile.c_str());
  fprintf(f,"%s\n",run->note.c_str());
  fprintf(f,"%s\n",run->suite.c_str());

  fprintf(f,"%d\n",(run->fiberRef[0]));
  fprintf(f,"%d\n",(run->fiberRef[1]));
  fprintf(f,"%d\n",(run->fiberRef[2]));
  fprintf(f,"%d\n",(run->fiberRef[3]));
  fprintf(f,"%d\n",(run->fiberRef[4]));
  fprintf(f,"%d\n",(run->fiberRef[5]));
  fprintf(f,"%d\n",(run->fiberRef[6]));
  fprintf(f,"%d\n",(run->fiberRef[7]));
  fprintf(f,"%d\n",(run->marocBoard[0]));
  fprintf(f,"%d\n",(run->marocBoard[1]));
  fprintf(f,"%d\n",(run->marocBoard[2]));
  fprintf(f,"%d\n",(run->marocBoard[3]));
  fprintf(f,"%d\n",(run->marocBoard[4]));
  fprintf(f,"%d\n",(run->marocBoard[5]));
  fprintf(f,"%d\n",(run->marocBoard[6]));
  fprintf(f,"%d\n",(run->marocBoard[7]));
  fprintf(f,"%d\n",(run->upstreamBoard));
  fprintf(f,"%f\n",(run->firstPath));
  fprintf(f,"%f\n",(run->secondPath));
  fprintf(f,"%f\n",(run->UpGEMz));
  fprintf(f,"%f\n",(run->DnGEMz));
  fprintf(f,"%f\n",(run->zAerogel));
  fprintf(f,"%lf\n",(run->geoCut));
  fprintf(f,"%lf\n",(run->radCut));
  fprintf(f,"%lf\n",(run->cutRadiusInRMS));
  fprintf(f,"%lf\n",(run->cutTimeInRMS));
  fprintf(f,"%lf\n",(run->cutRadiusOutRMS));
  fprintf(f,"%lf\n",(run->cutTimeOutRMS));
  fprintf(f,"%lf\n",(run->innerCorrectionX));
  fprintf(f,"%lf\n",(run->innerCorrectionY));
  fprintf(f,"%lf\n",(run->outerCorrectionX));
  fprintf(f,"%lf\n",(run->outerCorrectionY));
  fprintf(f,"%f\n",(run->UpGEMxRunOff));
  fprintf(f,"%f\n",(run->UpGEMyRunOff));
  fprintf(f,"%f\n",(run->DnGEMxRunOff));
  fprintf(f,"%f\n",(run->DnGEMyRunOff));
  fprintf(f,"%lf\n",(run->timeInMin));
  fprintf(f,"%lf\n",(run->timeInMax));
  fprintf(f,"%lf\n",(run->timeOuMin));
  fprintf(f,"%lf\n",(run->timeOuMax));
  fprintf(f,"%lf\n",(run->px474));
  fprintf(f,"%lf\n",(run->px519));
  fprintf(f,"%lf\n",(run->px537));
  fprintf(f,"%d\n",(run->beamChLogic));
  fprintf(f,"%d\n",(run->lookbackDAQ));
  fprintf(f,"%d\n",(run->aerogelRefractiveIndex));
  
  fprintf(f,"%lf\n",(run->MaxHitLength));
  fprintf(f,"%lf\n",(run->GlobalTimeOff));

  fclose(f);
  cout <<"Header recap wrote\n";
}




void writeHeaderExtended(THeader *run){
  FILE *f;
  string fname=Form("%s/output/header/run%04d_header.txt",run->suite.c_str(),run->runNum);
  f=fopen(fname.c_str(),"w");
  fprintf(f,"dRICH run number %d\n",run->runNum);
  fprintf(f,"GEM run number %d\n",run->runNumGEM);
  fprintf(f,"Date info: day %s, start time %s, end time %s\n",run->day.c_str(),run->startTime.c_str(),run->endTime.c_str());
  fprintf(f,"There was a beam of %s with energy %lf GeV\n",run->beam.c_str(),run->energyGeV);
  fprintf(f,"The aerogel mirror was in z=%f, the gas mirror was in z=%f\n",run->firstMirrorPosition,run->secondMirrorPosition);
  fprintf(f,"Data were acquired using the %ss\n",run->sensor.c_str());
  fprintf(f,"The trigger was %s\n",run->trigger.c_str());
  fprintf(f,"The reconstruction was performed using the %s setup file\n\n",run->setupFile.c_str());

  fprintf(f,"General run variables\n");
  fprintf(f,"Time coincidence window extremes In: %lf %lf Ou: %lf %lf\n",run->timeInMin,run->timeInMax,run->timeOuMin,run->timeOuMax);
  fprintf(f,"Upstream GEM on-beam centering offsets: %f %f\n",run->UpGEMxRunOff,run->UpGEMyRunOff);
  fprintf(f,"Downstream GEM on-beam centering offsets: %f %f\n",run->DnGEMxRunOff,run->DnGEMyRunOff); 

  fprintf(f,"Inner correction and cuts\n");
  fprintf(f,"x: %f\n",run->innerCorrectionX);
  fprintf(f,"y: %f\n",run->innerCorrectionY);
  fprintf(f,"Cut on radius: %lf\n",run->cutRadiusInRMS);
  fprintf(f,"Cut on time: %lf\n\n",run->cutTimeInRMS);
  fprintf(f,"Outer correction and cuts\n");
  fprintf(f,"x: %f\n",run->outerCorrectionX);
  fprintf(f,"y: %f\n",run->outerCorrectionY);
  fprintf(f,"Cut on radius: %lf\n",run->cutRadiusOutRMS);
  fprintf(f,"Cut on time: %lf\n",run->cutTimeOutRMS);

  fclose(f); 
}

