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

int FiberToPlace(int fiber){
  //PLACE LEGEND:
  //Sensore Nord = 0
  //Sensore Est = 1
  //Sensore Sud = 2
  //Sensore Ovest = 3
  //Fibra non riconosciuta, cout errore and return 5;
  int place = 5;
  if(fiber == 5 || fiber == 11){ //NORD
    place = 0;
  }else if( fiber == 6 || fiber == 9){ //EST
    place = 1;
  }else if( fiber == 4 || fiber == 8){ //SUD
    place = 2;
  }else if( fiber == 7 || fiber == 10){ //OVEST
    place = 3;
  }else{
    //cout <<Form("Errore nelle fibre!%d\n",fiber);
    //cin.get();
  }
  return place;
}


void MAPMTposition(int channel, int place, double *x, double *y){
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

}

/*void MPPCposition(int CHANNEL, int place,  double *x, double *y){
  int dX = (CHANNEL-1)/16;
  int dY = (CHANNEL-1)%16;
  if(place == 0 || place == 1){
    double iX=8.5*passoMPPC, iY=-7.5*passoMPPC;
    *x=iX+dX*passoMPPC;
    *y=iY+dY*passoMPPC;
  }
  if(place == 4 || place == 6){
    double iX=-7.5*passoMPPC, iY=-8.5*passoMPPC;
    *x=iX+dY*passoMPPC;
    *y=iY-dX*passoMPPC;
  }
  if(place == 5 || place == 7){
    double iX=-7.5*passoMPPC, iY=23.5*passoMPPC;
    *x=iX+dY*passoMPPC;
    *y=iY-dX*passoMPPC;
  }
}
*/
