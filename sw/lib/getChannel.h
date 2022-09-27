#include "definition.h"


void envVarCheck();
void getMaps();
int getMarocChip(int mCh);
int getMarocBoard(int fiber, THeader *run);
void upstreamMaroc(int fiber, THeader *run);
double timeCalibrationMAPMT(double time, int channel, int pmt);
double timeCalibrationMPPC(double time, int channel, int pmt);
int getMPPC_ch(int fiber, int mCh, int marocBoard, int chip, bool marocUpstream);
int getMAPMT_anode(int fiber, int mCh, int marocBoard, int chip, bool marocUpstream);
