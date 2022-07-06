#include <vector>

#include "dRICH.h"

using namespace std;

  //Default contructor and destructor;

dRICH::dRICH(){
  evtdRICH=-1;
  nedge=0;
  trigtime=0;  
  float x0=0;
  float y0=0;
  float cx0=0;
  float cy0=0;
  float x1=0;
  float y1=0;
  float cx1=0;
  float cy1=0; 
}

dRICH::~dRICH(){};

//Std method: write and read each member.
void dRICH::setEventdRICH(unsigned int Event){evtdRICH=Event;}
unsigned int dRICH::getEventdRICH(){return evtdRICH;}

void dRICH::setNedge(unsigned int Nedge){nedge=Nedge;}
unsigned int dRICH::getNedge() {return nedge;}


void dRICH::setData(int Slot, int Fiber, int ChM, int Pol, int Time){
  slot.push_back(Slot);
  fiber.push_back(Fiber);
  chM.push_back(ChM);
  pol.push_back(Pol);
  time.push_back(Time);
}
int dRICH::getSlot(int i){
  if((unsigned int)i < nedge)
    return slot[i];
  else
    return -1;
}

int dRICH::getFiber(int i){
  if((unsigned int)i < nedge)
    return fiber[i];
  else
    return -1;
}
int dRICH::getChM(int i){
  if((unsigned int)i < nedge)
    return chM[i];
  else
    return -1;
}
int dRICH::getPol(int i){
  if((unsigned int)i < nedge)
    return pol[i];
  else
    return -1;
}
int dRICH::getTime(int i){
  if((unsigned int)i < nedge)
    return time[i];
  else
    return -1;
}

void dRICH::setGEMs(float _x0, float _y0, float _cx0, float _cy0, float _x1, float _y1, float _cx1, float _cy1){
  x0=_x0;
  y0=_y0;
  cx0=_cx0;
  cy0=_cy0;
  x1=_x1;
  y1=_y1;
  cx1=_cx1;
  cy1=_cy1;
}

float dRICH::getX0(){return x0;}
float dRICH::getY0(){return y0;}
float dRICH::getCX0(){return cx0;}
float dRICH::getCY0(){return cy0;}
float dRICH::getX1(){return x1;}
float dRICH::getY1(){return y1;}
float dRICH::getCX1(){return cx1;}
float dRICH::getCY1(){return cy1;}
