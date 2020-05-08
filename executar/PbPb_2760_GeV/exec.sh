#!/bin/bash

#Number of the job in question
export NJOB=600
echo "NJOB= "$NJOB

#exporting eviroment variables necessary for the run
export FileDir="/sampa/fcanedo/IC/exechydro/PbPb_2760_GeV"
export MacroDir="/sampa/fcanedo/jewel-2.2.0"
export LD_LIBRARY_PATH="/sampa/fcanedo/lhapdf/lib"
export LHAPATH="/sampa/fcanedo/lhapdf/share/lhapdf/PDFsets"
export RIVET_ANALYSIS_PATH="/sampa/fcanedo/rivetanalises"
export HEPMC_FILE="/sampa/fcanedo/jewel-2.2.0/hepmc/PbPb_2760_GeV/cent0_10/exec."$NJOB".hepmc"
export YODA_GENERIC_FILE="/sampa/fcanedo/results/histograms/PbPb_2760_GeV/cent0_10/exec/MC_GENERIC/exec."$NJOB".yoda"
export YODA_SHAPE_FILE="/sampa/fcanedo/results/histograms/PbPb_2760_GeV/cent0_10/shape_exec/JET_SHAPE/shape.exec."$NJOB".yoda"
export YODA_ZHADRON_FILE="/sampa/fcanedo/results/histograms/PbPb_2760_GeV/cent0_10/z_hadron/exec."$NJOB".yoda"

chmod 777 $MacroDir/jewel-2.2.0-simple

source /cvmfs/alice.cern.ch/etc/login.sh

eval `alienv printenv VO_ALICE@GSL::v1.16-25`                                   
eval `alienv printenv VO_ALICE@boost::v1.59.0-21`                               
eval `alienv printenv VO_ALICE@HepMC::v2.06.09-13`                              
eval `alienv printenv VO_ALICE@fastjet::v3.2.1_1.024-alice3-1`                  
eval `alienv printenv VO_ALICE@YODA::v1.7.0-1`                                  
eval `alienv printenv VO_ALICE@Rivet::2.6.0-alice1-1`                           

cd $MacroDir/params/PbPb_2760_GeV/
cp exec.dat exec."$NJOB".dat

sed -i "s/\(NJOB\s\).*/\1"$NJOB"/g" exec."$NJOB".dat
sed -i "s/\(LOGFILE\s.*exec\.\).*/\1"$NJOB".log/g" exec."$NJOB".dat
sed -i "s/\(HEPMCFILE\s.*exec\.\).*/\1"$NJOB".hepmc/g" exec."$NJOB".dat
sed -i "s/\(MEDIUMPARAMS\s.*exec\.\).*/\1"$NJOB".dat/g" exec."$NJOB".dat

cd $MacroDir/medparams/PbPb_2760_GeV/
cp exec.dat exec."$NJOB".dat

sed -i "s/\(MEDFILE\s.*\/\)[0-9]*\(\.dat\)/\1"$MED_FILE_NUMBER"\2/g" exec."$NJOB".dat

cd $MacroDir
#$MacroDir/jewel-2.2.0-simple $MacroDir/params/PbPb_2760_GeV/exec."$NJOB".dat

cd $MacroDir/params/PbPb_2760_GeV/
rm exec."$NJOB".dat

cd $MacroDir/medparams/PbPb_2760_GeV/
rm exec."$NJOB".dat

export PSI=0.0

#rivet -a MC_GENERIC --ignore-beams -H $YODA_GENERIC_FILE $HEPMC_FILE
rivet  -a JET_SHAPE --ignore-beams -H $YODA_SHAPE_FILE $HEPMC_FILE
#rivet -a Z_HADRON --ignore-beams -H $YODA_ZHADRON_FILE $HEPMC_FILE
