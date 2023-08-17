#!/bin/bash

cmake --build .

sleep 2

for I in "274" "275" "276" "277" "278" "279" "280" "281"
do
  ./reco $I
  ./mon $I
done
