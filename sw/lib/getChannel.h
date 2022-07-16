#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

#include <TSystem.h>

#include "definition.h"

using namespace std;

int getMPAPMT_ch(int fiber, int mCh, int marocBoard, int chip);
void getMaps();
int getMarocBoard(int fiber, THeader *run);
int getMarocChip(int mCh);
string getMAPMT_ch(int fiber, int mChannel);
string getMAPMT_ch(int fiber, int mCh, int marocBoard, int chip);
string getMPPC_ch(int fiber, int mCh, int marocBoard, int chip);
