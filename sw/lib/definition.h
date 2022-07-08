#ifndef DEFINITION_H
#define DEFINITION_H

#include <string>
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

#endif
