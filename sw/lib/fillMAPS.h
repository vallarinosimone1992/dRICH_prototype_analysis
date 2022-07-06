#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>

//#include "headerStruct.h"

using namespace std;

struct header{
  int runNum;
  string day;
  string startTime;
  string endTime;
  string beam;
  int energyGeV;
  int expEvents;
  string sensor;
  float firstMirrorPosition;
  float secomdMirrorPosition;
  float temperature;
  int powerHV;
  string trigger;
  string runType;
  int runNumGEM;
  int pedestalGEM;
  string note;
};


void getRunNumbers(map<int,int> *map_GEM_run, map<int,double> *map_AeroMirrorPosition, map<int,double> *map_GasMirrorPosition);
void getMapMAPMT(map<string,int> *map_MAPMT1, map<string,int> *map_MAPMT2);
void readHeaders(int run, header *runHeader);
