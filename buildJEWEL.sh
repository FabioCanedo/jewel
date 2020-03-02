#!/bin/bash

source /cvmfs/alice.cern.ch/etc/login.sh

eval `alienv printenv VO_ALICE@lhapdf5::v5.9.1-2`                         

#alienv printenv VO_ALICE@lhapdf5::v5.9.1-alice3-1

make ; make clean
