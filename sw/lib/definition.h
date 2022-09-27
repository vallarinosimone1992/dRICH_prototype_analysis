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
  double geoCut=55.0;//55.0;
  double radCut=100.0;//100.0
  //Sigma of residui distribution for run 214
  double cutRadiusInRMS=4;// 1.93;
  double cutTimeInRMS=4;// 2.04;
  double cutRadiusOutRMS=4;// 2.56;
  double cutTimeOutRMS=4;// 2.59;
  
  double innerCorrectionX=0;
  double innerCorrectionY=0;
  double outerCorrectionX=0;
  double outerCorrectionY=0;
  float UpGEMxRunOff;
  float UpGEMyRunOff;
  float DnGEMxRunOff;
  float DnGEMyRunOff;
  double timeInMin;
  double timeInMax;
  double timeOuMin;
  double timeOuMax;
  double durMin=35;
  string outputDir;

  double px474=1.0;
  double px519=0.6;
  double px537=2.0;
  int beamChLogic;	
  int lookbackDAQ;

  double MaxHitLength=100;
  double GlobalTimeOff = 2400-lookbackDAQ;
};

static const bool SHOW_PROGRESS = false; 
static const int  CUT_MIN_DUR = 30; 
static const double CUT_NSIGMA_TIME = 3;    // nsigma for time selection
static const double CUT_FIT_COINC   = 5;    // width gate for coincidence fit (ns)

static const double RMS_ANGLE_GAS   = 1.4;  // min RMS in angle to exclude bad photons 
static const double RMS_TIME_GAS    = 3.0;  // min RMS iin time to exclude bad photons 
static const double RMS_ANGLE_AER   = 6.0;  // min RMS in angle to exclude bad photons 
static const double RMS_TIME_AER    = 3.0;  // min RMS iin time to exclude bad photons 

static const double GEM_CUT_X=20; //Maximum X for GEM
static const double GEM_CUT_Y=20; //Maximum Y for GEM
static const double GEM_CUT_R=.001; //Maximum theta for GEM

static const double mPi = 0.1396;
static const double mK = 0.49368; 
static const double mPr = 0.93827;

static const double nCO2 = 1.000410;
static const double nN2 = 1.000282;

//MAPMT
static const double xBinMAPMT[] = {-90,-81,-77.75, -74.71875, -71.6875, -68.65625, -65.625, -62.59375, -59.5625, -56.53125, -53.5, -50.46875, -47.4375, -44.40625, -41.375, -38.34375, -35.3125, -32.28125, -29.25, -24.25, -21.21875, -18.1875, -15.15625, -12.125, -9.09375, -6.0625, -3.03125, 0, 3.03125, 6.0625, 9.09375, 12.125, 15.15625, 18.1875, 21.21875, 24.25, 29.25, 32.28125, 35.3125, 38.34375, 41.375, 44.40625, 47.4375, 50.46875, 53.5, 56.53125, 59.5625, 62.59375, 65.625, 68.65625, 71.6875, 74.71875, 77.75,81,90};
static const double yBinMAPMT[] = {-90,-81,-77.75, -74.71875, -71.6875, -68.65625, -65.625, -62.59375, -59.5625, -56.53125, -53.5, -50.46875, -47.4375, -44.40625, -41.375, -38.34375, -35.3125, -32.28125, -29.25, -24.25, -21.21875, -18.1875, -15.15625, -12.125, -9.09375, -6.0625, -3.03125, 0, 3.03125, 6.0625, 9.09375, 12.125, 15.15625, 18.1875, 21.21875, 24.25, 29.25, 32.28125, 35.3125, 38.34375, 41.375, 44.40625, 47.4375, 50.46875, 53.5, 56.53125, 59.5625, 62.59375, 65.625, 68.65625, 71.6875, 74.71875, 77.75,81,90};
  //MPPC
static const double xBinMPPC[] = {-90,-85,-81.6,-78.4,-75.2,-72.0,-68.8,-65.6,-62.4,-59.2,-56.0,-52.8,-49.6,-46.4,-43.2,-40.0,-36.8,-33.6,-30.4,-25.6,-22.4,-19.2,-16.0,-12.8,-9.6,-6.4,-3.2,0.0,3.2,6.4,9.6,12.8,16.0,19.2,22.4,25.6,30.4,33.6,36.8,40.0,43.2,46.4,49.6,52.8,56.0,59.2,62.4,65.6,68.8,72.0,75.2,78.4,81.6,85,90};
static const double yBinMPPC[] = {-90,-85,-81.6,-78.4,-75.2,-72.0,-68.8,-65.6,-62.4,-59.2,-56.0,-52.8,-49.6,-46.4,-43.2,-40.0,-36.8,-33.6,-30.4,-25.6,-22.4,-19.2,-16.0,-12.8,-9.6,-6.4,-3.2,0.0,3.2,6.4,9.6,12.8,16.0,19.2,22.4,25.6,30.4,33.6,36.8,40.0,43.2,46.4,49.6,52.8,56.0,59.2,62.4,65.6,68.8,72.0,75.2,78.4,81.6,85,90};



#endif
