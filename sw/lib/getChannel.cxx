#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>
#include <string>

//#include <TLatex.h>
#include <TSystem.h>

#include "getChannel.h"
#include "fillMAPS.h"

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
  getRunNumbers(&m1,&m2,&m3);
  getMapMAPMT(&m4,&m5);
  getMapMPPC(&m6,&m7);
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


int getMAPMT_ch(int fiber, int mCh, int marocBoard, int chip, bool marocUpstream){
  if(m4.empty()){
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
  if(marocUpstream) nCh = m4.at(label);
  else  nCh = m5.at(label);
  if(nCh==-1){
    cout <<"[ERROR] Something wrong linking setup file and maps. Check it!\n";
    exit(EXIT_FAILURE);
  }
  return nCh;
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

