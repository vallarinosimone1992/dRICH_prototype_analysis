#define MAXDATA 10000
#include <iostream>
#include <stdio.h>
#include <string>
#include <map>
#include <iterator>
#include <vector>

#include <TObject.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TEllipse.h>
#include <TAxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TError.h>
#include <TMath.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TRandom3.h>
#include <TDirectory.h>

#include "../lib/fillMAPS.h"
#include "../lib/getChannel.h"
#include "../lib/MAPMTposition.h"
#include "../lib/handleFile.h"
#include "../lib/integrate.h"
#include "../lib/dRICH.h"


using namespace std;

map<int,int> m1;
map<int,int>::iterator it_m1;

map<int,double> m2;
map<int,double> m3;
map<int,double>::iterator it_m2;
map<int,double>::iterator it_m3;

map<string,int> m4;
map<string,int> m5;
map<string,int>::iterator it_m4;
map<string,int>::iterator it_m5;

const double xBin[] = {-81,-77.75, -74.71875, -71.6875, -68.65625, -65.625, -62.59375, -59.5625, -56.53125, -53.5, -50.46875, -47.4375, -44.40625, -41.375, -38.34375, -35.3125, -32.28125, -29.25, -24.25, -21.21875, -18.1875, -15.15625, -12.125, -9.09375, -6.0625, -3.03125, 0, 3.03125, 6.0625, 9.09375, 12.125, 15.15625, 18.1875, 21.21875, 24.25, 29.25, 32.28125, 35.3125, 38.34375, 41.375, 44.40625, 47.4375, 50.46875, 53.5, 56.53125, 59.5625, 62.59375, 65.625, 68.65625, 71.6875, 74.71875, 77.75,81};
const double yBin[] = {-81,-77.75, -74.71875, -71.6875, -68.65625, -65.625, -62.59375, -59.5625, -56.53125, -53.5, -50.46875, -47.4375, -44.40625, -41.375, -38.34375, -35.3125, -32.28125, -29.25, -24.25, -21.21875, -18.1875, -15.15625, -12.125, -9.09375, -6.0625, -3.03125, 0, 3.03125, 6.0625, 9.09375, 12.125, 15.15625, 18.1875, 21.21875, 24.25, 29.25, 32.28125, 35.3125, 38.34375, 41.375, 44.40625, 47.4375, 50.46875, 53.5, 56.53125, 59.5625, 62.59375, 65.625, 68.65625, 71.6875, 74.71875, 77.75,81};


int main(int argc, char *argv[]){
  TApplication theApp("App",&argc,argv);
  TCanvas *c = new TCanvas();
  c->Draw();
  TH2D *h = new TH2D("h","Rings; x[mm]; y[mm]",sizeof(xBin)/sizeof(*xBin)-1,xBin,sizeof(yBin)/sizeof(*yBin)-1,yBin);
  getRunNumbers(&m1, &m2, &m3);
  getMapMAPMT(&m4,&m5);

  TFile *fOut = new TFile("out.root","RECREATE");
  //TDirectory *dir = gDirectory();
  TTree *T = new TTree("dRICH","dRICH");
  TTreeIntegration(144,948,T);
  cout <<"Integration done\n";
  fOut->cd();
  T->Write();
  fOut->Close();
  return 0;
  
  
  //RANDOM FILL THE TH2D h
  TRandom3 rnd;
  rnd.SetSeed(01234567);
  for(int i = 0; i < 1e5; i++){
    int fiber = 4+rnd.Integer(8);
    int mCh = rnd.Integer(192);
    int channel=-1;
    string label = getMAPMT_ch(fiber, mCh);
    cout <<Form("Tentativo %d-esimo:\nFiber %d mCh %d\n",i+1,fiber,mCh);
    cout <<label <<endl <<endl;
    if(label.compare("N/A")==0) continue;
    if(fiber == 4 || fiber == 5 || fiber == 6 || fiber == 7){
      channel = m4.at(label);
    }else if(fiber == 8 || fiber == 9 || fiber == 10 || fiber == 11){
      channel = m5.at(label);
    }
    cout <<"channel: " <<channel <<endl;
    double x=0,y=0;
    int place = FiberToPlace(fiber);
    MAPMTposition(channel, place,&x, &y);
    cout <<Form("MAPMT ch position: %f,%f\n",x,y);
    h->Fill(x,y);
    //cin.get();
  }
  h->Draw("colz");
  c->Update();
  //theApp.Run();
  return 0;
}
