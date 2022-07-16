#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

//#include <TLatex.h>
#include <TSystem.h>

#include "getChannel.h"

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


void getMaps(){
  cout <<"Filling the static maps\n";
  getRunNumbers(&m1,&m2,&m3);
  getMapMAPMT(&m4,&m5);
  getMapMPPC(&m6,&m7);
}

bool upstreamMaroc(int fiber, THeader *run){
  for(int i = 0; i < 8; i++){
    if(fiber == run->fiberRef[i]){
      break;
    }
  }
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

/*
   int getMAPMT_ch(int fiber, int mCh, map<string,int> map_MAPMT1, map<string,int> map_MAPMT2){
   string label;
   int channel = -1;
   if(fiber == 4 || fiber == 5 || fiber == 6 || fiber == 7){
//MAROC 2 boards 
if(mChannel < 64){
label=Form("IN1_%02d",mChannel);
channel=map_MAPMT1.at(label.c_str());
}
if(mChannel>127){
label=Form("IN3_%02d",mChannel-128);
channel=map_MAPMT1.at(label.c_str());
}
}else if(fiber == 8 || fiber == 9 || fiber == 10 || fiber == 11){
//MAROC 3 boards
if(mChannel > 63 && mChannel < 128){
label=Form("IN1_%02d",mChannel-64);
channel=map_MAPMT2.at(label.c_str());
}
if(mChannel>127){
label=Form("IN3_%02d",mChannel-128);
channel=map_MAPMT2.at(label.c_str());
}
}else{
cout <<"Unknown board type! Press enter to continue";
cin.get();
}
return channel;
}
*/


int getMPAPMT_ch(int fiber, int mCh, int marocBoard, int chip){
  if(m4.empty())getMaps();
  string label = "N/A";
  int nCh=-1;
  if(chip==0) label=Form("IN1_%02d",mCh);
  else if(chip==1) label=Form("IN1_%02d",mCh-64);
  else if(chip==2) label=Form("IN3_%02d",mCh-128);
  else {
    cout <<Form("[ERROR] Wrong chip number, marocBoard 2, chip %d\n",chip);
    exit(EXIT_FAILURE);
  }
  if(marocUpstream) nCh = map4.at(label);
  else  nCh = map5.at(label);
  if(nCh==-1){
    cout <<"[ERROR] Something wrong linking setup file and maps. Check it!\n";
    exit(EXIT_FAILURE);
  }
  return nCh;
}

string getMAPMT_ch(int fiber, int mCh, int marocBoard, int chip){
  string label="N/A";
  if(marocBoard==2)
    if(chip==0)label=Form("IN1_%02d",mCh);
    else if(chip==2)label=Form("IN3_%02d",mCh-128);
    else{
      cout <<"[ERROR] something wrong in maroc channel setup!\n";
      cout <<fiber <<" " <<mCh <<" " <<marocBoard <<" " <<chip <<endl;
      exit(EXIT_FAILURE);
    }
  if(marocBoard==3){
    if(chip==1)label=Form("IN1_%02d",mCh-64);
    else if(chip==2)label=Form("IN3_%02d",mCh-128);
    else{
      cout <<"[ERROR] something wrong in maroc channel setup!\n";
      cout <<fiber <<" " <<mCh <<" " <<marocBoard <<" " <<chip <<endl;
      exit(EXIT_FAILURE);
    }
  }
  return label;
}

string getMPPC_ch(int fiber, int mCh, int marocBoard, int chip){
  string label="N/A";
  if(marocBoard==2){
    //if(mCh < 64)label = Form("IN1_%d",mCh);
    //else if(mCh > 127)label = Form("IN3_%d",mCh-128);
    if(chip==0)label = Form("IN1_%d",mCh);
    else if(chip==2)label = Form("IN3_%d",mCh-128);
    else{
      cout <<"[ERROR] something wrong in maroc channel setup!\n";
      cout <<fiber <<" " <<mCh <<" " <<marocBoard <<" " <<chip <<endl;
      exit(EXIT_FAILURE);
    }
  }
  if(marocBoard==3){
    //if(mCh > 63 && mCh < 128) label=Form("IN1_%d",mCh-64);
    //else if(mCh>127) label=Form("IN3_%d",mCh-128);
    if(chip==1) label=Form("IN1_%d",mCh-64);
    else if(chip==2) label=Form("IN3_%d",mCh-128);
    else{
      cout <<"[ERROR] something wrong in setup!\n";
      cout <<fiber <<" " <<mCh <<" " <<marocBoard <<" " <<chip <<endl;
      exit(EXIT_FAILURE);
    }
  }
  return label;
}


/*
   string label;
   int channel=300;
   if(fpga == 0 || fpga == 4 || fpga == 5){
   if(mChannel < 64){
   label=Form("IN1_%d",mChannel);
   channel=map_MPPC1.at(label.c_str());

   }
   if(mChannel>127){
   label=Form("IN3_%d",mChannel-128);
   channel=map_MPPC1.at(label.c_str());
   }
   }
   if(fpga == 1 || fpga == 6 || fpga == 7){
   if(mChannel < 64){
   label=Form("IN1_%d",mChannel);
   channel=map_MPPC2.at(label.c_str());
   }
   if(mChannel>127){
   label=Form("IN3_%d",mChannel-128);
   channel=map_MPPC2.at(label.c_str());
   }
   }
   return channel;

*/

