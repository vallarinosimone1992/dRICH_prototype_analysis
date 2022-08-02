#ifndef UTILITY_H
#define UTILITY_H
#include <iostream>
#include "definition.h"

void printProgress(double progress);
void printEnd();

void printUsageReco();
void printUsageAna();

void checkFileExistance(THeader *run);

void writeHeader(THeader *run);
#endif
