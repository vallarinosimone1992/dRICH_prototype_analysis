#!/bin/bash

cmake --build .

sleep 2

if [ $1 -gt 0 ]
then
	./reco 127		
	./reco 130		
	./reco 131		
	./reco 132		
	./reco 133
fi

./mon $2 127 130 131 132 133
