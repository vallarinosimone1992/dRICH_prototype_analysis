#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

#include <TSystem.h>

#include "definition.h"

using namespace std;

void envVarCheck();
void getMaps();
int getMarocChip(int mCh);
int getMarocBoard(int fiber, THeader *run);
void upstreamMaroc(int fiber, THeader *run);
int getMPPC_ch(int fiber, int mCh, int marocBoard, int chip, bool marocUpstream);
int getMAPMT_ch(int fiber, int mCh, int marocBoard, int chip, bool marocUpstream);
