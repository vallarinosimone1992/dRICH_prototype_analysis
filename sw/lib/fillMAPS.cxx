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
  file=fopen("../map/MAPMT.map","r");
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


void readHeaders(int run, THeader *runHeader){
  FILE *file;
  file=fopen("../../../DATA/header/logbook.tsv","r");
  char line0[10000];
  char line[10000];
  int tRunNum, tEnergyGeV, tExpEvents, tPowerHV, tRunNumGEM, tPedestalGEM;
  char tDay[200], tStartTime[200],tEndTime[200],tBeam[200], tSensor[200], tTrigger[200],tRunType[200],tNote[2000];
  float tFirstMirrorPosition, tSecondMirrorPosition, tTemperature;
  auto n =fgets(line0,10000,file);
  while(fgets(line,10000,file)!=NULL){
    sscanf(line,"%d %s %s %s %s %d %d %s %f %f %f %d %s %s %d %d %s",&tRunNum,tDay,tStartTime,tEndTime,tBeam,&tEnergyGeV,&tExpEvents,tSensor,&tFirstMirrorPosition,&tSecondMirrorPosition,&tTemperature,&tPowerHV,tTrigger,tRunType,&tRunNumGEM,&tPedestalGEM,tNote);
    if(run == tRunNum){
      cout <<"The header includes the following info\n";
      cout <<line0;
      cout <<line;
      //cout <<Form("%d %s %s %s %s %d %d %s %f %f %f %d %s %s %d %d %s\n",tRunNum,&tDay[0],&tStartTime[0],&tEndTime[0],&tBeam[0],tEnergyGeV,tExpEvents,&tSensor[0],tFirstMirrorPosition,tSecondMirrorPosition,tTemperature,tPowerHV,&tTrigger[0],&tRunType[0],tRunNumGEM,tPedestalGEM,&tNote[0]);
      runHeader->runNum = tRunNum;
      runHeader->day=tDay;
      runHeader->startTime=tStartTime;
      runHeader->endTime=tEndTime;
      runHeader->beam=&tBeam[0];
      runHeader->energyGeV=tEnergyGeV;
      runHeader->expEvents=tExpEvents;
      runHeader->sensor=tSensor;
      runHeader->firstMirrorPosition=tFirstMirrorPosition;
      runHeader->secomdMirrorPosition=tSecondMirrorPosition;
      runHeader->temperature=tTemperature;
      runHeader->powerHV=tPowerHV;
      runHeader->trigger=tTrigger;
      runHeader->runType=tRunType;
      runHeader->runNumGEM=tRunNumGEM;
      runHeader->pedestalGEM=tPedestalGEM;
      runHeader->note=tNote;

      break;
    }
  }
}
