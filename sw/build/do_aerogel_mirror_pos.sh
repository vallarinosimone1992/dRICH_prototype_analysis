#!/bin/bash

cmake --build .

sleep 2

for I in "107" "108" "109" "110" "111" "112" "113" "114" "115" "116" "117" "118"
do
  ./reco $I
  ./mon $I
done
