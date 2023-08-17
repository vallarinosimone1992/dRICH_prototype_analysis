#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

//#include <TLatex.h>
#include <TSystem.h>
#include <TH1D.h>

#include "getChannel.h"
#include "fillMAPS.h"
#include "correction.h"

using namespace std;

static map<int,int> m1;
static map<int,int>::iterator it_m1;

static map<int,double> m2;
static map<int,double> m3;
static map<int,double>::iterator it_m2;
static map<int,double>::iterator it_m3;

static map<string,int> m4;
static map<string,int> m5;
static map<string,int>::iterator it_m4;
static map<string,int>::iterator it_m5;

static map<string,int> m6;
static map<string,int> m7;
static map<string,int>::iterator it_m6;
static map<string,int>::iterator it_m7;

static map<int,double> m8;
static map<int,double> m9;
static map<int,double>::iterator it_m8;
static map<int,double>::iterator it_m9;


static vector<double> qeLambda;
static vector<double> qeEfficiency;

void getMaps(){
  getRunNumbers(&m1,&m2,&m3);
  getMapMAPMT(&m4,&m5);
  getMapMPPC(&m6,&m7);
  getTimeCalibrationDataMAPMT(&m8);
  getTimeCalibrationDataMPPC(&m9);
  getQuantumEfficiencyValues(&qeLambda,&qeEfficiency);
}

void getQuantumEfficiencyValues(vector<double> *qeLambda, vector<double> *qeEfficiency){
  FILE *file;
  file=fopen("../map/QE_value.dat","r");
  char line[400];
  double lambda, eff;
  while(fgets(line,400,file)!=NULL){
    sscanf(line,"%lf %lf",&lambda,&eff);
    qeLambda->push_back(lambda);
    qeEfficiency->push_back(eff/100);
  }
}

bool applyQuantumEfficiency(double energy){
  //if(qeLambda.size() < 1)getQuantumEfficiencyValues();
  double cmp = 1240./energy/1e6;
  int label=-1;
  bool out=false;
  if(rnd.Rndm() > 0.8)return false;
  for(int i = 0; i < qeLambda.size()-1; i++){
    if(qeLambda[i] > cmp && qeLambda[i+1] <= cmp){
      label = i;
      break;
    }
  }
  if(label < 0) return false;
  double random = rnd.Rndm();
  double cmpEff=(qeEfficiency[label]+qeEfficiency[label+1])/2;
  if(random >cmpEff)out=false;
  else out=true;
  return out;
}



void upstreamMaroc(int fiber, THeader *run){
  bool answer=false;
  for(int i = 0; i < 8; i+=2){
    if(fiber == run->fiberRef[i]){
      answer=true;
      break;
    }
  }
  run->upstreamBoard=answer;
}

int getMarocBoard(int fiber, THeader *run){
  int board = -1;
  for(int i = 0; i < 8; i++){
    if(fiber == run->fiberRef[i]){
      board = run->marocBoard[i];
      break;
    }
  } 
  return board;
}

int getMarocChip(int mCh){
  int ret = -1;
  if(mCh < 64) ret = 0;
  if(mCh > 63 && mCh < 128) ret = 1;
  if(mCh > 127) ret = 2;
  return ret;
}


double timeCalibrationMAPMT(double time, int channel, int pmt){
  if(m8.empty())getMaps();
  int place = pmt*256+channel;
  double correction = m8.at(place);
  //double newTime = time + correction;
  double newTime;
  if(correction!=0) newTime = time - correction+390; //Nicola calibration
  else newTime = time - correction; //Nicola calibration
                                    //if(channel<10)cout <<"Place: "<<pmt<<place<<" Correcion:  "<<time <<" + "<<correction <<" = " <<newTime <<endl;
                                    //if(channel<10)printf("Place (%2d, %3d)  %4d corr %7.2f %7.2f = %7.2f \n",pmt,channel,place,time,correction,newTime);
                                    //cin.get();
  return newTime;
}

double timeCalibrationMPPC(double time, int channel, int pmt){
  if(m9.empty())getMaps();
  int place = pmt*256+channel;
  double correction = m9.at(place);
  double newTime = time + correction;
  return newTime;
}


int getMAPMT_anode(int fiber, int mCh, int marocBoard, int chip, bool marocUpstream){
  if(m4.empty()){
    //cout <<"[WARNING] filling again the static maps\n";
    getMaps();
  }
  string label = "N/A";
  int anode=-1;
  if(marocBoard==2 && chip==0) label=Form("IN1_%02d",mCh);
  else if(marocBoard==3 && chip==1) label=Form("IN1_%02d",mCh-64);
  else if(chip==2) label=Form("IN3_%02d",mCh-128);
  else if(marocBoard==3 && chip==0) anode=-1; 
  else {
    cout <<Form("[ERROR] Wrong chip number, marocBoard 2, chip %d\n",chip);
    exit(EXIT_FAILURE);
  }
  if(marocBoard!=3 || chip!=0){
    if(marocUpstream) anode = m4.at(label);
    else  anode = m5.at(label);
    if(anode==-1){
      cout <<"[ERROR] Something wrong linking setup file and maps. Check it!\n";
      exit(EXIT_FAILURE);
    }
  }
  return anode;
}

int getMPPC_ch(int fiber, int mCh, int marocBoard, int chip, bool marocUpstream){
  if(m6.empty()){
    //cout <<"[WARNING] filling again the static maps\n";
    getMaps();
  }
  string label = "N/A";
  int nCh=-1;
  if(chip==0) label=Form("IN1_%02d",mCh);
  else if(chip==1) label=Form("IN1_%02d",mCh-64);
  else if(chip==2) label=Form("IN3_%02d",mCh-128);
  else {
    cout <<Form("[ERROR] Wrong chip number, marocBoard 2, chip %d\n",chip);
    exit(EXIT_FAILURE);
  }
  if(marocUpstream) nCh=m6.at(label);
  else nCh=m7.at(label);
  if(nCh==-1){
    cout <<"[ERROR] Something wrong linking setup file and maps. Check it!\n";
    exit(EXIT_FAILURE);
  }
  return nCh; 
}
