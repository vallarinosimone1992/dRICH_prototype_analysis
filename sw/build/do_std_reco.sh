#!/bin/bash

echo "If you want to do the reconstruction again, modify the script!"

sleep 2

cmake --build .

sleep 2

#if [ $# -eq 0 ]
  #then
    #./reco 130
    #./reco 131
    #./reco 132
    #./reco 133
#fi

./mon std 130 131 132 133
