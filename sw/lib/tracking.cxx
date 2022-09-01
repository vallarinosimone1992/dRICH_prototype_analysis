#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>

#include "tracking.h"
#include "definition.h"

//GEM constant values provided by Evaristo.
const double offset=51.2;
const double pitch_GEM=0.4;

double GEM_getBeamlineOffset(TH1D *h){
  gErrorIgnoreLevel=kWarning;
  TF1 *f = new TF1("f","gaus(0)",-100,100);
  f->SetParameters(1200,0,4);
  f->SetParLimits(0,1,5000);
  h->Fit("f","Q","N",-40,40);
  h->Draw();
  return f->GetParameter(1);
}

void GEM_relative(float *x0, float *y0, float *x1, float *y1){
  *x0 = *x0 * pitch_GEM-offset;
  *y0 = *y0 * pitch_GEM-offset;
  *x1 = *x1 * pitch_GEM-offset;
  *y1 = *y1 * pitch_GEM-offset;
}

void GEM_position(THeader *run, float *x0, float *y0, float *x1, float *y1){
  *x0=*x0-run->UpGEMxRunOff; 
  *y0=*y0-run->UpGEMyRunOff; 
  *x1=*x1-run->DnGEMxRunOff; 
  *y1=*y1-run->DnGEMyRunOff; 
  //cout <<"Position "<<*x0 <<" " <<*y0 <<" " <<*x1 <<" " <<*y1 <<endl;
}


void AERO_computing(THeader *run, float *xA, float *yA, float *mx, float *my, float x0, float y0, float x1, float y1){
  *mx = (x1-x0)/(run->DnGEMz-run->UpGEMz);
  *my = (y1-y0)/(run->DnGEMz-run->UpGEMz);
  float qx = x0;
  float qy = y0;
  *xA= *mx * run->zAerogel + qx;
  *yA= *my * run->zAerogel + qy;
  //cout <<"Aerogel: " <<*xA <<" " <<*yA <<" " <<*mx <<" " <<*my <<endl;
}

