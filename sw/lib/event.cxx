#include <vector>

#include "event.h"

using namespace std;

  //Default contructor and destructor;

event::event(){
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
  float gx0=0;
  float gy0=0;
  float gx1=0;
  float gy1=0;
  float thetaX=0;
  float thetaY=0;
  float aeroX=0;
  float aeroY=0;
}

event::~event(){};

//Std method: write and read each member.
void event::setEventdRICH(unsigned int Event){evtdRICH=Event;}
unsigned int event::getEventdRICH(){return evtdRICH;}

void event::setNedge(unsigned int Nedge){nedge=Nedge;}
unsigned int event::getNedge() {return nedge;}


void event::setData(int Slot, int Fiber, int ChM, int Pol, int Time){
  slot.push_back(Slot);
  fiber.push_back(Fiber);
  chM.push_back(ChM);
  pol.push_back(Pol);
  time.push_back(Time);
}
int event::getSlot(int i){
  if((unsigned int)i < nedge)
    return slot[i];
  else
    return -1;
}

int event::getFiber(int i){
  if((unsigned int)i < nedge)
    return fiber[i];
  else
    return -1;
}
int event::getChM(int i){
  if((unsigned int)i < nedge)
    return chM[i];
  else
    return -1;
}
int event::getPol(int i){
  if((unsigned int)i < nedge)
    return pol[i];
  else
    return -1;
}
int event::getTime(int i){
  if((unsigned int)i < nedge)
    return time[i];
  else
    return -1;
}

void event::setGEMs(float _x0, float _y0, float _cx0, float _cy0, float _x1, float _y1, float _cx1, float _cy1){
  x0=_x0;
  y0=_y0;
  cx0=_cx0;
  cy0=_cy0;
  x1=_x1;
  y1=_y1;
  cx1=_cx1;
  cy1=_cy1;
}

float event::getX0(){return x0;}
float event::getY0(){return y0;}
float event::getCX0(){return cx0;}
float event::getCY0(){return cy0;}
float event::getX1(){return x1;}
float event::getY1(){return y1;}
float event::getCX1(){return cx1;}
float event::getCY1(){return cy1;}


void event::setGEMsDerived(){}

float event::getGX0(){return gx0;}
float event::getGY0(){return gy0;}
float event::getGX1(){return gx1;}
float event::getGY1(){return gy1;}
