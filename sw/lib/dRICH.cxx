#include <vector>

#include "dRICH.h"

using namespace std;

  //Default contructor and destructor;

dRICH::dRICH(){
  evtdRICH=-1;
  nedge=0;
  trigtime=0;  
   float x0;
  float y0;
  float cx0;
  float cy0;
  float x1;
  float y1;
  float cx1;
  float cy1;
 
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
