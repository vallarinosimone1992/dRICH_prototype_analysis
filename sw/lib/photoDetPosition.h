#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

#include <TSystem.h>

using namespace std;

int FiberToPlace(int fiber);
void MAPMTposition(int channel, int place, double *x, double *y);
void MPPCposition(int channel, int place, double *x, double *y);
void SIPMposition(int channel, int place, double *x, double *y);
