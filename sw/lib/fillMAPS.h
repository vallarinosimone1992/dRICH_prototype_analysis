#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>

#include "definition.h"

using namespace std;


void getRunNumbers(map<int,int> *map_GEM_run, map<int,double> *map_AeroMirrorPosition, map<int,double> *map_GasMirrorPosition);
void getMapMAPMT(map<string,int> *map_MAPMT1, map<string,int> *map_MAPMT2);
void getMapMPPC(map<string,int> *map_MPPC1, map<string,int> *map_MPPC2);
void readHeaders(int run, THeader *runHeader);
