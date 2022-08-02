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
  if(fgets(line,MAXLINEA,f))sscanf(line,"%d",&(run->energyGeV));
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
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->timeMin));
  if(fgets(line,MAXLINEA,f))sscanf(line,"%lf",&(run->timeMax));
  fclose(f);
}


void writeHeaderShort(THeader *run){
  ofstream fout;
  string fname=Form("%s/output/header/run%04d_header_data.txt",run->suite.c_str(),run->runNum);
  fout.open(fname.c_str());

  fout <<run->runNum<<endl;
  fout <<run->day.c_str()<<endl;
  fout <<run->startTime.c_str()<<endl;
  fout <<run->endTime.c_str()<<endl;
  fout <<run->beam<<endl;
  fout <<run->energyGeV<<endl;
  fout <<run->expEvents<<endl;
  fout <<run->sensor<<endl;
  fout <<run->firstMirrorPosition<<endl;
  fout <<run->secondMirrorPosition<<endl;
  fout <<run->temperature<<endl;
  fout <<run->powerHV<<endl;
  fout <<run->trigger<<endl;
  fout <<run->runType<<endl;
  fout <<run->runNumGEM<<endl;
  fout <<run->pedestalGEM<<endl;
  fout <<run->setupFile<<endl;
  fout <<run->note<<endl;
  fout <<run->suite<<endl;
  fout <<run->fiberRef[0]<<endl;
  fout <<run->fiberRef[1]<<endl;
  fout <<run->fiberRef[2]<<endl;
  fout <<run->fiberRef[3]<<endl;
  fout <<run->fiberRef[4]<<endl;
  fout <<run->fiberRef[5]<<endl;
  fout <<run->fiberRef[6]<<endl;
  fout <<run->fiberRef[7]<<endl;
  fout <<run->marocBoard[0]<<endl;
  fout <<run->marocBoard[1]<<endl;
  fout <<run->marocBoard[2]<<endl;
  fout <<run->marocBoard[3]<<endl;
  fout <<run->marocBoard[4]<<endl;
  fout <<run->marocBoard[5]<<endl;
  fout <<run->marocBoard[6]<<endl;
  fout <<run->marocBoard[7]<<endl;
  fout <<run->upstreamBoard<<endl;
  fout <<run->firstPath<<endl;
  fout <<run->secondPath<<endl;
  fout <<run->UpGEMz<<endl;
  fout <<run->DnGEMz<<endl;
  fout <<run->zAerogel<<endl;
  fout <<run->geoCut<<endl;
  fout <<run->cutRadiusInRMS<<endl;
  fout <<run->cutTimeInRMS<<endl;
  fout <<run->cutRadiusOutRMS<<endl;
  fout <<run->cutTimeOutRMS<<endl;
  fout <<run->innerCorrectionX<<endl;
  fout <<run->innerCorrectionY<<endl;
  fout <<run->outerCorrectionX<<endl;
  fout <<run->outerCorrectionY<<endl;
  fout <<run->UpGEMxRunOff<<endl;
  fout <<run->UpGEMyRunOff<<endl;
  fout <<run->DnGEMxRunOff<<endl;
  fout <<run->DnGEMyRunOff<<endl;
  fout <<run->timeMin<<endl;
  fout <<run->timeMax<<endl;
  fout.close();
  cout <<"Header recap wrote\n";
}



void writeHeaderExtended(THeader *run){
  FILE *f;
  string fname=Form("%s/output/header/run%04d_header.txt",run->suite.c_str(),run->runNum);
  f=fopen(fname.c_str(),"w");
  fprintf(f,"dRICH run number %d\n",run->runNum);
  fprintf(f,"GEM run number %d\n",run->runNumGEM);
  fprintf(f,"Date info: day %s, start time %s, end time %s\n",run->day.c_str(),run->startTime.c_str(),run->endTime.c_str());
  fprintf(f,"There was a beam of %s with energy %d GeV\n",run->beam.c_str(),run->energyGeV);
  fprintf(f,"The aerogel mirror was in z=%f, the gas mirror was in z=%f\n",run->firstMirrorPosition,run->secondMirrorPosition);
  fprintf(f,"Data were acquired using the %ss\n",run->sensor.c_str());
  fprintf(f,"The trigger was %s\n",run->trigger.c_str());
  fprintf(f,"The reconstruction was performed using the %s setup file\n\n",run->setupFile.c_str());

  fprintf(f,"General run variables\n");
  fprintf(f,"Time coincidence window extremes: %lf %lf\n",run->timeMin,run->timeMax);
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

