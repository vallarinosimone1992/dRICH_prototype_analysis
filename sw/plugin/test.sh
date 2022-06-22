#!/bin/bash

clear

rm exe_test.exe

g++ -I `root-config --incdir --cflags` -Wall -O2 ../lib/* main_test.C -o exe_test.exe `root-config --libs`

./exe_test.exe
