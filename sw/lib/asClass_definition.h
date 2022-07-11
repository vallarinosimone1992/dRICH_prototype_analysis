#ifndef DEFINITION_H
#define DEFINITION_H

#include <TObject.h>
#include <TCollection.h>
#include <string>
using namespace std;

class  THeader : public TObject{
  private:
    int runNum;
    string day;
    string startTime;
    string endTime;
    string beam;
    int beamEnergy;
    int expEvents;
    string sensor;
    float firstMirrorPosition;
    float secomdMirrorPosition;
    float temperature;
    int powerHV;
    string trigger;
    string runType;
    int runNumGEM;
    int pedestalGEM;
    string note;
  public:
    THeader();
    THeader(char *line[]);
    int RunNum();
    string Day();
    string StartTime();
    string EndTime();
    string Beam();
    int BeamEnergy();
    int ExpEvents();
    string Sensor();
    float 


};

#endif
