#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>

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

