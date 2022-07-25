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
  float secondMirrorPosition;
  float temperature;
  int powerHV;
  string trigger;
  string runType;
  int runNumGEM;
  int pedestalGEM;
  string setupFile;
  string note;
  string suite;
  int fiberRef[8]={-1,-1,-1,-1,-1,-1,-1,-1};
  int marocBoard[8]={-1,-1,-1,-1,-1,-1,-1,-1};
  bool upstreamBoard;
  float firstPath;
  float secondPath;
  float UpGEMz;
  float DnGEMz;
  float zAerogel;
  double geoCut=55.0;
  double innerCorrectionX=0;
  double innerCorrectionY=0;
  double outerCorrectionX=0;
  double outerCorrectionY=0;
  float UpGEMxRunOff;
  float UpGEMyRunOff;
  float DnGEMxRunOff;
  float DnGEMyRunOff;
  double timeMin;
  double timeMax;
};

#endif
