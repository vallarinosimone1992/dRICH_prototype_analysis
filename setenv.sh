#!/bin/bash

if [[ -z "${DRICH_SUITE}" ]]; then
  echo "You have to set the environment variables DRICH_SUITE as the upstream directory"
  exit
fi

echo "The env variable DRICH_SUITE is $DRICH_SUITE"

mkdir -p $DRICH_SUITE/DATA/dRICH_DATA
mkdir -p $DRICH_SUITE/DATA/GEM_DATA
mkdir -p $DRICH_SUITE/DATA/header
mkdir -p $DRICH_SUITE/processed_data/firstStepData
mkdir -p $DRICH_SUITE/processed_data/integrated_dRICH_GEM_data
mkdir -p $DRICH_SUITE/output/header
mkdir -p $DRICH_SUITE/output/plot

echo "Environment set completed"


