#include "definition.h"

double mmTomRad(double r, double inPath, double zMir);
void convertToRadiant(THeader *run);
void applyFit(TH1D *h, TF1 *f,string fname, bool out);
void computeRMS(THeader *run);
void computeCutSingleParticle(THeader *run);
void newSingleParticle(THeader *run);
void singleParticle(THeader *run);
