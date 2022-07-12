#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

#include <TSystem.h>

using namespace std;

int FiberToPhDet(int fiber, int cmp[8]);
int FiberToPlace(int fiber);
int FiberToMAPMT(int fiber);
int FiberToMPPC(int fiber);
void MAPMTposition(int channel, int place, double *x, double *y, double *r);
void MPPCposition(int channel, int place, double *x, double *y, double *r);
void MAPMTposition(int channel, int place, double *x, double *y);
void SIPMposition(int channel, int place, double *x, double *y);
