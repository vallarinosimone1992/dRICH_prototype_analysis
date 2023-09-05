#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>

#include "tracking.h"
#include "definition.h"


//GEM constant values provided by Evaristo.
const double offset=51.2;
const double pitch_GEM=0.4;

double GEM_getBeamlineOffset(TH1D *h, int pri){
  gErrorIgnoreLevel=kWarning;
  TF1 *f = new TF1("f","gaus(0)",-100,100);
  f->SetParameters(1200,h->GetBinCenter(h->GetMaximumBin()),5);
  if(pri==1)printf("Fit Y1 %7.2f \n",h->GetBinCenter(h->GetMaximumBin()));
  //f->SetParLimits(2,1,10);
  h->Fit("f","","",h->GetBinCenter(h->GetMaximumBin())-40,h->GetBinCenter(h->GetMaximumBin())+40);
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  h->Draw();
  //if(pri)c1->Print("gy1.pdf");
  return f->GetParameter(1);
  //return h->GetBinCenter(h->GetMaximumBin());
}

void GEM_relative(float *x0, float *y0, float *x1, float *y1){
  *x0 = *x0 * pitch_GEM-offset;
  *x1 = *x1 * pitch_GEM-offset;
  if(ANALYSIS_2022==true){
    *y0 = *y0 * pitch_GEM-offset;
    *y1 = *y1 * pitch_GEM-offset;
  }else{
    *y0 = -(*y0 * pitch_GEM-offset);
    *y1 = -(*y1 * pitch_GEM-offset);
  }
}

void GEM_position(THeader *run, float *x0, float *y0, float *x1, float *y1){
  *x0=(*x0-run->UpGEMxRunOff); 
  *y0=(*y0-run->UpGEMyRunOff); 
  *x1=(*x1-run->DnGEMxRunOff); 
  *y1=(*y1-run->DnGEMyRunOff); 
  //cout <<"Position "<<*x0 <<" " <<*y0 <<" " <<*x1 <<" " <<*y1 <<endl;
}


void AERO_computing(THeader *run, float *xA, float *yA, float *mx, float *my, float x0, float y0, float x1, float y1){
  *mx = (x1-x0)/(run->DnGEMz-run->UpGEMz);
  *my = (y1-y0)/(run->DnGEMz-run->UpGEMz);
  float qx = x0;
  float qy = y0;
  *xA= *mx * run->zAerogel + qx;
  *yA= *my * run->zAerogel + qy;
  *mx = atan2((x1-x0),(run->DnGEMz-run->UpGEMz));
  *my = atan2((y1-y0),(run->DnGEMz-run->UpGEMz));
}

