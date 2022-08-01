#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>

#include "definition.h"

static TH1D* hTime;

void inizializePlot(){
  hTime = new TH1D("hTime","Time distribution for all hit",1000,0,1000);
}

void fillTime(double time){
  hTime->Fill(time);
}

void plotAll(){
  TCanvas *c1 = new TCanvas("c1","c1",1600,900);
  c1->Draw();
  hTime->Draw();
  c1->Update();
}
