#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

#include <TSystem.h>

#include "photoDetPosition.h"

using namespace std;

double minRMAPMT  = 29.25;
double frameMAPMT = (51.8-48.5)/2.;  // PMT frame structure
double halfMAPMT  = 48.5/2.;         // active area
double firstMAPMT = (48.5-42)/2.;    // first pixel
double passoMAPMT = 3.0;             // inner pixels

double min1MPPC=30.4;
double max1MPPC=81.6;
double min2MPPC=25.6;
double passoMPPC = 3.2;

//double Losanga = 54.25;
double Losanga = 0;



//---------------------------------------------------
void MAPMTposition(int channel, int place, double *x, double *y, double *r){
  //---------------------------------------------------
  double row = (channel-1)/16;
  double col = (channel-1)%16;
  double step = passoMAPMT;
  //if(row==0 || row==15) step = (firstMAPMT+passoMAPMT)/2.;
  //if(col==0 || col==15) step = (firstMAPMT+passoMAPMT)/2.;
  double step0 = (firstMAPMT+passoMAPMT)/2;

  if(place == 0){//NORD
    //double iX=-halfMAPMT+0.5*firstMAPMT, iY=minRMAPMT+frameMAPMT+0.5*firstMAPMT;
    double iX=-halfMAPMT-Losanga+0.5*firstMAPMT, iY=minRMAPMT+0.5*firstMAPMT;
    if(row == 0) *x = iX;
    else if(row == 15) *x=iX+2*step0+(row-2)*step;
    else *x = iX+step0+(row-1)*step;
    //*x=iX+row*step;
     if(col == 0) *y = iY;
    else if(col == 15) *y=iY+2*step0+(col-2)*step;
    else *y = iY+step0+(col-1)*step;
    //*y=iY+col*step;
  }
  if(place == 1){//EST
    //double iX=minRMAPMT+frameMAPMT+0.5*firstMAPMT, iY=halfMAPMT-0.5*firstMAPMT;
    //*x=iX+col*step;
    //*y=iY-row*step;
    double iX=minRMAPMT+0.5*firstMAPMT, iY=halfMAPMT+Losanga-0.5*firstMAPMT;
    if(col==0)*x=iX;
    else if(col==15)*x=iX+2*step0+(col-2)*step;
    else *x = iX+step0+(col-1)*step;
    if(row==0)*y=iY;
    else if(row==15)*y=iY-2*step0-(row-2)*step;
    else *y=iY-step0-(row-1)*step;
  }
  if(place == 2){//SUD
    //double iX=halfMAPMT-0.5*firstMAPMT, iY=-minRMAPMT-frameMAPMT-0.5*firstMAPMT;
    //*x=iX-row*step;
    //*y=iY-col*step;
    double iX=halfMAPMT+Losanga-0.5*firstMAPMT, iY=-minRMAPMT-0.5*firstMAPMT;
    if(row==0)*x=iX;
    else if(row==15)*x=iX-2*step0-(row-2)*step;
    else *x=iX-step0-(row-1)*step;
    if(col==0) *y=iY;
    else if(col==15) *y=iY-2*step0-(col-2)*step;
    else *y=iY-step0-(col-1)*step;
  }
  if(place == 3){//WEST
    //double iX=-minRMAPMT-frameMAPMT-0.5*firstMAPMT, iY=-halfMAPMT+0.5*firstMAPMT;
    //*x=iX-col*step;
    //*y=iY+row*step;
    double iX=-minRMAPMT-0.5*firstMAPMT, iY=-halfMAPMT-Losanga+0.5*firstMAPMT;
    if(col==0)*x=iX;
    else if(col==15)*x=iX-2*step0-(col-2)*step;
    else *x=iX-step0-(col-1)*step;
    if(row==0) *y=iY;
    else if(row==15)*y=iY+2*step0+(row-2)*step;
    else *y=iY+step0+(row-1)*step;
  }
  *r = sqrt((*x)*(*x)+(*y)*(*y));
}




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

void MPPCposition(int CHANNEL, int place,  double *x, double *y, double *r){
  int ratio = (CHANNEL-1)/16;
  int remainder = (CHANNEL-1)%16;
  //Place tab:
  //0 -> North. Pixel 1 in bottom right
  //1 -> Est
  //2 -> South
  if(place == 0){
    //double iX=-7.5*passoMPPC, iY=23.5*passoMPPC;
    double iX=-min2MPPC+.5*passoMPPC, iY=max1MPPC-0.5*passoMPPC;
    *x=iX+remainder*passoMPPC;
    *y=iY-ratio*passoMPPC;
  }
  if(place == 1){
    //double iX=8.5*passoMPPC, iY=-7.5*passoMPPC;
    double iX=min1MPPC+0.5*passoMPPC, iY=-min2MPPC+0.5*passoMPPC;
    *x=iX+ratio*passoMPPC;
    *y=iY+remainder*passoMPPC;
  }
  if(place == 2){
    //double iX=-7.5*passoMPPC, iY=-8.5*passoMPPC;
    double iX=-min2MPPC+0.5*passoMPPC, iY=-min1MPPC-0.5*passoMPPC;
    *x=iX+remainder*passoMPPC;
    *y=iY-ratio*passoMPPC;
    //    cout <<"Channel " <<CHANNEL <<" " <<ratio <<" " <<remainder <<"; X: " <<iX <<" " <<*x <<"; Y: " <<iY <<" " <<*y <<endl;
    //    cin.get();
  }
  *r = sqrt((*x)*(*x)+(*y)*(*y));
}

