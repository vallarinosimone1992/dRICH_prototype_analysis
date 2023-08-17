#!/bin/bash

cmake --build .

sleep 2

if [ $# -eq 0 ]
  then
    ./reco 140
    ./reco 142
fi

./mon 50GeV 140 142 
