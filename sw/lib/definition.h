#ifndef DEFINITION_H
#define DEFINITION_H

#include <string>
#include <TString.h>
using namespace std;

struct THeader{
  int runNum;
  string day;
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
  //Sigma of residui distribution for run 214
  double cutRadiusInRMS=10;// 1.93;
  double cutTimeInRMS=10;// 2.04;
  double cutRadiusOutRMS=10;// 2.56;
  double cutTimeOutRMS=10;// 2.59;
  
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
