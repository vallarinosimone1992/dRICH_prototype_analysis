#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

#include <TSystem.h>

#include "photoDetPosition.h"

using namespace std;

double min1=29.25;
double min2=24.25;
double passoMAPMT = 3.03125;
double passoMPPC = 3.2;

int FiberToPhDet(int fiber, int cmp[8]){
  //PLACE LEGEND:
  //North = 0, Est = 1, South = 2, West = 3
  //Looking the photo-detection plane
  int place = 5;
  if(fiber == cmp[0] || fiber == cmp[1]){ //NORTH
    place = 0;
  }else if(fiber == cmp[2] || fiber == cmp[3]){ //EST
    place = 1;
  }else if(fiber == cmp[4] || fiber == cmp[5]){ //SOUTH
    place = 2;
  }else if(fiber == cmp[6] || fiber == cmp[7]){ //WEST
    place = 3;
  }else{
    //cout <<Form("Errore nelle fibre!%d\n",fiber);
    //cin.get();
  }
  return place;
}


void MAPMTposition(int channel, int place, double *x, double *y, double *r){
  double dX = (channel-1)/16;
  double dY = (channel-1)%16;
  if(place == 0){//NORD
    double iX=-(min2-0.5*passoMAPMT), iY=min1+.5*passoMAPMT;
    *x=-(iX+dX*passoMAPMT);
    *y=iY+dY*passoMAPMT;
  }
  if(place == 1){//EST
    double iX=min1+0.5*passoMAPMT, iY=min2-.5*passoMAPMT;
    *x=iX+dY*passoMAPMT;
    *y=-(iY-dX*passoMAPMT);
  }
  if(place == 2){//SUD
    double iX=min2-.5*passoMAPMT, iY=-(min1+.5*passoMAPMT);
    *x=-(iX-dX*passoMAPMT);
    *y=iY-dY*passoMAPMT;
  }
  if(place == 3){//WEST
    double iX=-(min1+.5*passoMAPMT), iY=(-min2+.5*passoMAPMT);
    *x=iX-dY*passoMAPMT;
    *y=-(iY+dX*passoMAPMT);
  }
  *r = sqrt((*x)*(*x)+(*y)*(*y));
}

void MPPCposition(int CHANNEL, int place,  double *x, double *y, double *r){
  int dX = (CHANNEL-1)/16;
  int dY = (CHANNEL-1)%16;
  //Place tab:
  //0 -> North
  //1 -> Est
  //2 -> South
  if(place == 0){
    double iX=-7.5*passoMPPC, iY=23.5*passoMPPC;
    *x=iX+dY*passoMPPC;
    *y=iY-dX*passoMPPC;
  }
  if(place == 1){
    double iX=8.5*passoMPPC, iY=-7.5*passoMPPC;
    *x=iX+dX*passoMPPC;
    *y=iY+dY*passoMPPC;
  }
  if(place == 2){
    double iX=-7.5*passoMPPC, iY=-8.5*passoMPPC;
    *x=iX+dY*passoMPPC;
    *y=iY-dX*passoMPPC;
  }
  *r = sqrt((*x)*(*x)+(*y)*(*y));
}

