#!/bin/bash

#Number of the job in question
export NJOB=163
echo "NJOB= "$NJOB

#exporting eviroment variables necessary for the run
#export FileDir="/sampa/fcanedo/IC/vusphydro/PbPb_5020_GeV"
export MacroDir="/sampa/fcanedo/jewel-2.2.0"
export LD_LIBRARY_PATH="/sampa/fcanedo/lhapdf/lib"
export LHAPATH="/sampa/fcanedo/lhapdf/share/lhapdf/PDFsets"
export RIVET_ANALYSIS_PATH="/sampa/fcanedo/rivetanalises"
export HEPMC_FILE="/sampa/fcanedo/jewel-2.2.0/hepmc/pp_2760_GeV/pp."$NJOB".hepmc"
export YODA_GENERIC_FILE="/sampa/fcanedo/results/histograms/pp_2760_GeV/pp/MC_GENERIC/pp."$NJOB".yoda"
export YODA_SHAPE_FILE="/sampa/fcanedo/results/histograms/pp_2760_GeV/shape_pp/JET_SHAPE/shape.pp."$NJOB".yoda"
export YODA_ZHADRON_FILE="/sampa/fcanedo/results/histograms/pp_2760_GeV/z_hadron/pp."$NJOB".yoda"

chmod 777 $MacroDir/jewel-2.2.0-simple

source /cvmfs/alice.cern.ch/etc/login.sh

eval `alienv printenv VO_ALICE@GSL::v1.16-25`                                   
eval `alienv printenv VO_ALICE@boost::v1.59.0-21`                               
eval `alienv printenv VO_ALICE@HepMC::v2.06.09-13`                              
eval `alienv printenv VO_ALICE@fastjet::v3.2.1_1.024-alice3-1`                  
eval `alienv printenv VO_ALICE@YODA::v1.7.0-1`                                  
eval `alienv printenv VO_ALICE@Rivet::2.6.0-alice1-1`                           

#export MED_FILE_NUMBER=`cat $FileDir/first_filelist.txt | head -n $NJOB | tail -n 1`
#echo "MED FILE NUMBER= "$MED_FILE_NUMBER

#cd $FileDir
#sed -i "13s/\(BEGIN=\).*/\1$NJOB/g" reader.sh
#sed -i "14s/\(END=\).*/\1$NJOB/g" reader.sh
#./reader.sh

#if [ -e "$FileDir"/"$MED_FILE_NUMBER".dat.xz ]; then
#xz -d "$FileDir"/"$MED_FILE_NUMBER".dat.xz
#echo "Medfile exists"
#else
#echo "Medfile does not exist"
#fi

cd $MacroDir/params/pp_5020_GeV/
cp pp.dat pp."$NJOB".dat

sed -i "s/\(NJOB\s\).*/\1"$NJOB"/g" pp."$NJOB".dat
sed -i "s/\(LOGFILE\s.*pp\.\).*/\1"$NJOB".log/g" pp."$NJOB".dat
sed -i "s/\(HEPMCFILE\s.*pp\.\).*/\1"$NJOB".hepmc/g" pp."$NJOB".dat
sed -i "s/\(MEDIUMPARAMS\s.*pp\.\).*/\1"$NJOB".dat/g" pp."$NJOB".dat

cd $MacroDir/medparams/pp_5020_GeV/
cp pp.dat pp."$NJOB".dat

sed -i "s/\(MEDFILE\s.*\/\)[0-9]*\(\.dat\)/\1"$MED_FILE_NUMBER"\2/g" pp."$NJOB".dat

cd $MacroDir
#$MacroDir/jewel-2.2.0-vac $MacroDir/params/pp_5020_GeV/pp."$NJOB".dat

cd $MacroDir/params/pp_5020_GeV/
rm pp."$NJOB".dat

cd $MacroDir/medparams/pp_5020_GeV/
rm pp."$NJOB".dat

#if [ ! -e "$FileDir"/"$MED_FILE_NUMBER".dat.xz ]; then
#echo ""
#xz "$FileDir"/"$MED_FILE_NUMBER".dat
#fi



#let line=($NJOB%100+1)
#export PSI=`cat jakiPsi | head -n $NJOB | tail -n 1`
#echo $PSI

#rivet -a MC_GENERIC --ignore-beams -H $YODA_GENERIC_FILE $HEPMC_FILE
#rivet  -a JET_SHAPE --ignore-beams -H $YODA_SHAPE_FILE $HEPMC_FILE
rivet -a Z_HADRON --ignore-beams -H $YODA_ZHADRON_FILE $HEPMC_FILE
