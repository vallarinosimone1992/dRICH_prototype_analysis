#define MAXDATA 10000
#include <iostream>
#include <stdio.h>
#include <map>
#include <iterator>
#include <vector>

#include <TTree.h>
#include <TMath.h>
#include <TSystem.h>


using namespace std;

void step1(TTree *t, TTree *tGEM, TTree *tout);
void TTreeIntegration(int runDRICH, int runGEM, TTree *tout);
void TTreeIntegration(int runDRICH, int runGEM);

/*int MAPMT_channel_three(int fpga, int mChannel){
  string label;
  int channel=300;
  if(fpga == 4 || fpga == 5 || fpga == 6 || fpga == 7){
    if(mChannel > 63 && mChannel < 128){
      label=Form("IN1_%02d",mChannel-64);
      channel=map_MAPMT1.at(label.c_str());
    }
    if(mChannel>127){
      label=Form("IN3_%02d",mChannel-128);
      channel=map_MAPMT1.at(label.c_str());
    }
  }
  if(fpga == 8 || fpga == 9 || fpga == 10 || fpga == 11){
    if(mChannel > 63 && mChannel < 128){
      label=Form("IN1_%02d",mChannel-64);
      channel=map_MAPMT2.at(label.c_str());
    }
    if(mChannel>127){
      label=Form("IN3_%02d",mChannel-128);
      channel=map_MAPMT2.at(label.c_str());
    }
  }
  return channel;
}

int MAPMT_channel_two(int fpga, int mChannel){
  string label;
  int channel=300;
  if(fpga == 4 || fpga == 5 || fpga == 6 || fpga == 7){
    if(mChannel < 64){
      label=Form("IN1_%02d",mChannel);
      channel=map_MAPMT1.at(label.c_str());
    }
    if(mChannel>127){
      label=Form("IN3_%02d",mChannel-128);
      channel=map_MAPMT1.at(label.c_str());
    }
  }
  if(fpga == 8 || fpga == 9 || fpga == 10 || fpga == 11){
    if(mChannel < 64){
      label=Form("IN1_%02d",mChannel);
      channel=map_MAPMT2.at(label.c_str());
    }
    if(mChannel>127){
      label=Form("IN3_%02d",mChannel-128);
      channel=map_MAPMT2.at(label.c_str());
    }
  }
  return channel;
}


int MPPC_channel_two(int fpga, int mChannel){
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
}
*/
