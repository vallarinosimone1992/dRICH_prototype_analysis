#ifndef DEFINITION_H
#define DEFINITION_H

#include <string>
#include <TString.h>
#include <TRandom3.h>
using namespace std;

static const bool debug=false;
static const bool APPLY_QUANTUM_EFFICIENCY=true;
static const bool APPLY_GEM_CUT=false;
static const bool APPLY_PIXELATION=true;
static const bool APPLY_SIMULATION_TRACKING_ERROR=false;
static const bool SWAP_UPSTREAM_DOWNSTREAM_GEM=true;
static const bool ANALYSIS_2022=true;
static TRandom3 rnd;

struct THeader{
  int runNum;
  string day;
  string startTime;
  string endTime;
  string beam;
  double energyGeV;
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
  int mergedRunCross=-1;
  int mergedRunCorner=-1;
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
  double aerogelRefractiveIndex=-1;
  
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
  double durMin=30;
  string outputDir;

  double px474=1.0;
  double px519=0.6;
  double px537=2.0;
  int beamChLogic;	
  int lookbackDAQ;

  double MaxHitLength=100;
  double GlobalTimeOff = 2400-lookbackDAQ;

  int beamLUND=211;
};

static const bool DATA_2021 = false;

static const bool SHOW_PROGRESS = false; 
static const int  CUT_MIN_DUR = 30; 
static const double CUT_NSIGMA_TIME = 3;    // nsigma for time selection
static const double CUT_FIT_COINC   = 5;    // width gate for coincidence fit (ns)

static const double RMS_ANGLE_GAS   = 1.4;  // min RMS in angle to exclude bad photons 
static const double RMS_TIME_GAS    = 3.0;  // min RMS iin time to exclude bad photons 
static const double RMS_ANGLE_AER   = 6.0;  // min RMS in angle to exclude bad photons 
static const double RMS_TIME_AER    = 3.0;  // min RMS iin time to exclude bad photons 

static const double GEM_CUT_X=10; //Maximum X for GEM [mm]
static const double GEM_CUT_Y=10; //Maximum Y for GEM [mm]
static const double GEM_CUT_R=.0005; //Maximum theta for GEM [Rad]

static const double mPi = 0.1396;
static const double mK = 0.49368; 
static const double mPr = 0.93827;

static const double nCO2 = 1.000410;
static const double nN2 = 1.000282;

static const double SIMULATION_DETECTOR_Z=-16; //Std position: -16

static const double SIMULATION_AERO_MIRROR_RADIUS=700;
static const double SIMULATION_AERO_MIRROR_CENTER=-380;
static const double SIMULATION_GAS_MIRROR_RADIUS=2420;
static const double SIMULATION_GAS_MIRROR_CENTER=-380;

static const double SIMULATION_FILTER_WIDTH=0;
static const double SIMULATION_AEROGEL_WIDTH=0;
static const double SIMULATION_AEROGEL_EXIT_Z=61;
static const double SIMULATION_MIRROR_Z=317;
static const double SIMULATION_MIRROR_HOLE_RADIUS=57.5;
//static const double SIMULATION_AERO_PATH=sqrt(pow(SIMULATION_MIRROR_HOLE_RADIUS,2)+pow(SIMULATION_MIRROR_Z,2))+SIMULATION_AEROGEL_EXIT_Z+SIMULATION_FILTER_WIDTH+SIMULATION_AEROGEL_WIDTH;
static const double SIMULATION_AERO_PATH=368;//Std position: 368
static const double SIMULATION_GAS_PATH=1130;//Std position 1130

static const double SIMULATION_GAS_PHOTON_EXPECTED_TIME=2.7;
static const double SIMULATION_AEROGEL_PHOTON_EXPECTED_TIME=8.5;

static const double UPSTREAM_GEM_TIME_SIM = 0;
static const double UPSTREAM_GEM_Z_SIM = -150;
static const double DOWNSTREAM_GEM_TIME_SIM = 6.7;
static const double DOWNSTREAM_GEM_Z_SIM = +1850;


//MAPMT
static const double xBinMAPMT[] = {-90,-77.75,-74.5,-71.5,-68.5,-65.5,-62.5,-59.5,-56.5,-53.5,-50.5,-47.5,-44.5,-41.5,-38.5,-35.5,-32.5,-29.25,-24.25,-21,-18,-15,-12,-9,-6,-3,0,3,6,9,12,15,18,21,24.25,29.25,32.5,35.5,38.5,41.5,44.5,47.5,50.5,53.5,56.5,59.5,62.5,65.5,68.5,71.5,74.5,77.75,90};
static const double yBinMAPMT[] = {-90,-77.75,-74.5,-71.5,-68.5,-65.5,-62.5,-59.5,-56.5,-53.5,-50.5,-47.5,-44.5,-41.5,-38.5,-35.5,-32.5,-29.25,-24.25,-21,-18,-15,-12,-9,-6,-3,0,3,6,9,12,15,18,21,24.25,29.25,32.5,35.5,38.5,41.5,44.5,47.5,50.5,53.5,56.5,59.5,62.5,65.5,68.5,71.5,74.5,77.75,90};



//static const double xBinMAPMT[] = {-90,-81,-77.75, -74.71875, -71.6875, -68.65625, -65.625, -62.59375, -59.5625, -56.53125, -53.5, -50.46875, -47.4375, -44.40625, -41.375, -38.34375, -35.3125, -32.28125, -29.25, -24.25, -21.21875, -18.1875, -15.15625, -12.125, -9.09375, -6.0625, -3.03125, 0, 3.03125, 6.0625, 9.09375, 12.125, 15.15625, 18.1875, 21.21875, 24.25, 29.25, 32.28125, 35.3125, 38.34375, 41.375, 44.40625, 47.4375, 50.46875, 53.5, 56.53125, 59.5625, 62.59375, 65.625, 68.65625, 71.6875, 74.71875, 77.75,81,90};
//static const double yBinMAPMT[] = {-90,-81,-77.75, -74.71875, -71.6875, -68.65625, -65.625, -62.59375, -59.5625, -56.53125, -53.5, -50.46875, -47.4375, -44.40625, -41.375, -38.34375, -35.3125, -32.28125, -29.25, -24.25, -21.21875, -18.1875, -15.15625, -12.125, -9.09375, -6.0625, -3.03125, 0, 3.03125, 6.0625, 9.09375, 12.125, 15.15625, 18.1875, 21.21875, 24.25, 29.25, 32.28125, 35.3125, 38.34375, 41.375, 44.40625, 47.4375, 50.46875, 53.5, 56.53125, 59.5625, 62.59375, 65.625, 68.65625, 71.6875, 74.71875, 77.75,81,90};
  //MPPC
static const double xBinMPPC[] = {-90,-85,-81.6,-78.4,-75.2,-72.0,-68.8,-65.6,-62.4,-59.2,-56.0,-52.8,-49.6,-46.4,-43.2,-40.0,-36.8,-33.6,-30.4,-25.6,-22.4,-19.2,-16.0,-12.8,-9.6,-6.4,-3.2,0.0,3.2,6.4,9.6,12.8,16.0,19.2,22.4,25.6,30.4,33.6,36.8,40.0,43.2,46.4,49.6,52.8,56.0,59.2,62.4,65.6,68.8,72.0,75.2,78.4,81.6,85,90};
static const double yBinMPPC[] = {-90,-85,-81.6,-78.4,-75.2,-72.0,-68.8,-65.6,-62.4,-59.2,-56.0,-52.8,-49.6,-46.4,-43.2,-40.0,-36.8,-33.6,-30.4,-25.6,-22.4,-19.2,-16.0,-12.8,-9.6,-6.4,-3.2,0.0,3.2,6.4,9.6,12.8,16.0,19.2,22.4,25.6,30.4,33.6,36.8,40.0,43.2,46.4,49.6,52.8,56.0,59.2,62.4,65.6,68.8,72.0,75.2,78.4,81.6,85,90};



#endif
