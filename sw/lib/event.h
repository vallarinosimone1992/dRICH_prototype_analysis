#include <iostream>
#include <vector>

using namespace std;

class event{
  private:
  const float x0_inOffset=51.2; //[mm]
  const float y0_inOffset=51.2;
  const float x1_inOffset=51.2;
  const float y1_inOffset=51.2;
  const float pitchGEM=0.4;



  public:
  unsigned int evtdRICH;
  unsigned int nedge;
  vector<int> slot;
  vector<int> fiber;
  vector<int> chM;
  vector<int> pol;
  vector<int> time;
  double trigtime;
  int row;
  int inst;
  int evtGEM;
  float x0;
  float y0;
  float cx0;
  float cy0;
  float x1;
  float y1;
  float cx1;
  float cy1;
  float gx0;
  float gy0;
  float gx1;
  float gy1;
  float thetaX;
  float thetaY;
  float aeroX;
  float aeroY;



  //Default contructor and destructor;
  event();
  virtual ~event();

  //Std method: write and read each member.
  void setEventdRICH(unsigned int Event);
  unsigned int getEventdRICH();

  void setNedge(unsigned int Nedge);
  unsigned int getNedge();

  void setData(int Slot, int Fiber, int ChM, int Pol, int Time);
  int getSlot(int i);
  int getFiber(int i);
  int getChM(int i);
  int getPol(int i);
  int getTime(int i);

  void setGEMs(float _x0, float _y0,float _cx0, float _cy0,float _x1, float _y1,float _cx1, float _cy1);
  float getX0();
  float getY0();
  float getCX0();
  float getCY0();
  float getX1();
  float getY1();
  float getCX1();
  float getCY1();

  void setGEMsDerived();
  float getGX0();
  float getGY0();
  float getGX1();
  float getGY1();
  float getThetaX();
  float getThetaY();
  float getAeroX();
  float getAeroY();
};

