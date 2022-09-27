#include <TH1D.h>
#include "definition.h"


double GEM_getBeamlineOffset(TH1D *h, int pri);
void GEM_relative(float *x0, float *y0, float *x1, float *y1);
void GEM_position(THeader *run, float *x0, float *y0, float *x1, float *y1);
void AERO_computing(THeader *run, float *xA, float *yA, float *mx, float *my, float x0, float y0, float x1, float y1);

