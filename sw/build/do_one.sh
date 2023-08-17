#!/bin/bash

cmake --build .

sleep 2

if [ $1 -lt 1999 ]
then
./reco $1
else
./sim $1
fi

./mon $2 $1
