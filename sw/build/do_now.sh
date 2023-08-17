#!/bin/bash


cmake --build .

sleep 2

for I in {285..299}
do
    echo ./reco $I
done

./mon PS_10GeV+ 285 286
./mon PS_10GeV- 287 288
./mon PS_8GeV+ 289 290
./mon PS_8GeV- 291 292
./mon PS_6GeV+ 293 294
./mon PS_6GeV- 295 296
./mon PS_4GeV+ 297
./mon PS_4GeV- 298
