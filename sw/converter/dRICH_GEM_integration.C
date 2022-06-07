#include <iostream>
#include <stdio.h>
#include <string>
#include <map>
#include <iterator>
#include <vector>

#include "../lib/fillMAPS.h"


using namespace std;

map<int,int> map_GEM_rNumber;
map<int,double> map_AeroMirror_Position;
map<int,double> map_GasMirror_Position;

map<int,int>::iterator it_GEM_rNumber;
map<int,double>::iterator it_AeroMirror_Position;
map<int,double>::iterator it_GasMirror_Position;

int main(int argc, char *argv[]){
  getRunNumbers(&map_GEM_rNumber,&map_AeroMirror_Position,&map_GasMirror_Position);
}
