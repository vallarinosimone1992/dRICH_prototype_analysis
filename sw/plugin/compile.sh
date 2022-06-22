#!/bin/bash

clear

g++ -I `root-config --incdir --cflags` -Wall -O2 ../lib/* dRICH_GEM_integration.C -o exe_dRICH_GEM_integration.exe `root-config --libs`

