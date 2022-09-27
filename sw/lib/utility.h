#ifndef UTILITY_H
#define UTILITY_H
#include <iostream>
#include "definition.h"

void printProgress(double progress);
void printEnd();

int reference(int i, int j);
void printUsageReco();
void printUsageAna();
void printUsageMon();

void checkFileExistance(THeader *run);

void writeHeader(THeader *run);
#endif
