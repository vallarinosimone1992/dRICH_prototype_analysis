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

string getMAPMT_ch(int fiber, int mChannel){
  string label="N/A";
  if(fiber == 4 || fiber == 5 || fiber == 6 || fiber == 7){
    //MAROC 2 boards 
    if(mChannel < 64){
      label=Form("IN1_%02d",mChannel);
    }
    if(mChannel>127){
      label=Form("IN3_%02d",mChannel-128);
    }
  }else if(fiber == 8 || fiber == 9 || fiber == 10 || fiber == 11){
    //MAROC 3 boards
    if(mChannel > 63 && mChannel < 128){
      label=Form("IN1_%02d",mChannel-64);
    }
    if(mChannel>127){
      label=Form("IN3_%02d",mChannel-128);
    }
  }else{
    cout <<"Unknown board type! Press enter to continue";
    cin.get();
  }
  return label;
}
