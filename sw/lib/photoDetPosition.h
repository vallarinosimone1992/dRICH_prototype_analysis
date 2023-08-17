#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

#include <TSystem.h>

using namespace std;

int FiberToPhDet(int fiber, int cmp[8]);
void MAPMTposition(int channel, int place, double *x, double *y, double *r);
void MPPCposition(int channel, int place, double *x, double *y, double *r);
bool simulationPixel(double x, double y, int *pmt, double *px, double *py, double *r);
