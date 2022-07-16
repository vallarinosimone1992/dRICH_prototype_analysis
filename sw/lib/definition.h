#ifndef DEFINITION_H
#define DEFINITION_H

#include <string>
#include <TString.h>
using namespace std;

struct THeader{
  int runNum;
  TString day;
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
  int fiberRef[8]={-1,-1,-1,-1,-1,-1,-1,-1};
  int marocBoard[8]={-1,-1,-1,-1,-1,-1,-1,-1};
  float firstPath;
  float secondPath;
};

#endif
