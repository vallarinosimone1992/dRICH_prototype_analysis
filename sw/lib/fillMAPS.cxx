#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>

#include <TSystem.h>

#include "fillMAPS.h"

void getRunNumbers(map<int,int> *map_GEM_run, map<int,double> *map_AeroMirrorPosition, map<int,double> *map_GasMirrorPosition){
  FILE *file;
  file=fopen("../map/run_number.map","r");
  char line[2000];
  int key;
  int gemRun;
  float gasPos, aeroPos;
  while(fgets(line,2000,file)!=NULL){
    sscanf(line,"%d %d %f %f",&key,&gemRun, &aeroPos, &gasPos);
    map_GEM_run->insert(make_pair(key,gemRun));
    map_AeroMirrorPosition->insert(make_pair(key,aeroPos));
    map_GasMirrorPosition->insert(make_pair(key,gasPos));
  }
  fclose(file);
}



void getMapMAPMT(map<string,int> *map_MAPMT1, map<string,int> *map_MAPMT2){
  FILE *file;
  file=fopen("../map/MAPMT_mapping_data.map","r");
  char line[2000];
  char key[2000];
  int value;
  while(fgets(line,2000,file)!=NULL){
    sscanf(line,"%d %s",&value,key);
    if( value <= 128)map_MAPMT1->insert(make_pair(key,value));
    if( value > 128)map_MAPMT2->insert(make_pair(key,value));
  }
  fclose(file);
}

void getMapMPPC(map<string,int> *map_MPPC1, map<string,int> *map_MPPC2){
  FILE *file;
  file=fopen("../map/MPPC_mapping_data.dat","r");
  if(file == NULL) cout <<"File not opened\n";
  char line[2000];
  char key[2000];
  int value;
  while(fgets(line,2000,file)!=NULL){
    sscanf(line,"%d %s",&value,key);
    if(value<=128)map_MPPC1->insert(make_pair(key,value));
    if(value>128)map_MPPC2->insert(make_pair(key,value));
  }
  fclose(file);
}

void getTimeCalibrationDataMAPMT(map<int,double> *map_time_MAPMT){
  FILE *file;
  file=fopen("../map/MAPMT_time_calibration.dat","r");
  if(file == NULL){
    cout <<"MAPMT time calibration map not found\n";
    exit(EXIT_FAILURE);
  }
  char line[2000];
  int key;
  double value;
  while(fgets(line,2000,file)!=NULL){
    sscanf(line,"%d %lf",&key,&value);
    map_time_MAPMT->insert(make_pair(key,value));
  }
  fclose(file);
}

void getTimeCalibrationDataMPPC(map<int,double> *map_time_MPPC){
  FILE *file;
  file=fopen("../map/MPPC_time_calibration.dat","r");
  if(file == NULL){ 
    cout <<"MPPC time calibration map not found\n";
    exit(EXIT_FAILURE);
  }
  char line[2000];
  int key;
  double value;
  while(fgets(line,2000,file)!=NULL){
    sscanf(line,"%d %lf",&key,&value);
    map_time_MPPC->insert(make_pair(key,value));
  }
  fclose(file);
}

void readHeaders(int run, THeader *runHeader){
  const char  *tmp = getenv("DRICH_SUITE");
  string env_var(tmp ? tmp : "");
  if(env_var.empty()){
    cerr <<"[ERROR] No such variable found! You should define the variable DRICH_SUITE!" <<endl;
    exit(EXIT_FAILURE);
  }
  runHeader->suite=env_var;
  FILE *file;
  string headerName=Form("%s/DATA/logbook/logbook.tsv",runHeader->suite.c_str());
  if(gSystem->AccessPathName(headerName.c_str())){
    cout <<Form("[ERROR] Header file %s not found\n",headerName.c_str());
    exit(EXIT_FAILURE);
  }
  file=fopen(headerName.c_str(),"r");
  char line0[10000];
  char line[10000];
  int tRunNum=-1, tEnergyGeV, tExpEvents, tPowerHV, tRunNumGEM, tPedestalGEM;
  char tDay[200], tStartTime[200],tEndTime[200],tBeam[200], tSensor[200], tTrigger[200],tRunType[200], tSetupFile[2000],tNote[2000];
  float tFirstMirrorPosition, tSecondMirrorPosition, tTemperature;
  auto n =fgets(line0,10000,file);
  bool headerFound=false;
  while(fgets(line,10000,file)!=NULL){
    sscanf(line,"%d %s %s %s %s %d %d %s %f %f %f %d %s %s %d %d %s %s",&tRunNum,tDay,tStartTime,tEndTime,tBeam,&tEnergyGeV,&tExpEvents,tSensor,&tFirstMirrorPosition,&tSecondMirrorPosition,&tTemperature,&tPowerHV,tTrigger,tRunType,&tRunNumGEM,&tPedestalGEM,tSetupFile,tNote);
    if(run == tRunNum){
      runHeader->runNum = tRunNum;
      runHeader->day=tDay;
      runHeader->startTime=tStartTime;
      runHeader->endTime=tEndTime;
      runHeader->beam=&tBeam[0];
      runHeader->energyGeV=tEnergyGeV;
      runHeader->expEvents=tExpEvents;
      runHeader->sensor=tSensor;
      runHeader->firstMirrorPosition=tFirstMirrorPosition;
      runHeader->secondMirrorPosition=tSecondMirrorPosition;
      runHeader->temperature=tTemperature;
      runHeader->powerHV=tPowerHV;
      runHeader->trigger=tTrigger;
      runHeader->runType=tRunType;
      runHeader->runNumGEM=tRunNumGEM;
      runHeader->pedestalGEM=tPedestalGEM;
      runHeader->setupFile=tSetupFile;
      runHeader->note=tNote;

      headerFound=true;
      break;
    }
  }
  fclose(file);

  if(headerFound==false){
    cout <<"[ERROR] Run number not found in the logbook\n";
    exit(EXIT_FAILURE);
  }

  file=fopen(Form("%s/dRICH_prototype_analysis/sw/map/%s.map",runHeader->suite.c_str(),runHeader->setupFile.c_str()),"r");
  for(int i = 0; i < 8; i++){
    if(fgets(line,10000,file)!=NULL){
      int tmp1, tmp2;
      auto prz = sscanf(line,"%d %d",&tmp1,&tmp2);
      runHeader->fiberRef[i]=tmp1;
      runHeader->marocBoard[i]=tmp2;
    }
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->firstPath=tmp1;
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->secondPath=tmp1;
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->UpGEMz=tmp1;
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->DnGEMz=tmp1;
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->zAerogel=tmp1;
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->geoCut=tmp1;
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->innerCorrectionX=tmp1;
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->innerCorrectionY=tmp1;
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->outerCorrectionX=tmp1;
  }
  if(fgets(line,10000,file)!=NULL){
    float tmp1;
    auto prz = sscanf(line,"%f",&tmp1);
    runHeader->outerCorrectionY=tmp1;
  }



  cout <<"Header and setup file read\n";
}
