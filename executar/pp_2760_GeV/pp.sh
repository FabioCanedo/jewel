#!/bin/bash

#Usage:
#$1=NJOB
#$2=NEVENT
#$3=energy
#$4=RUN


NEVENT=$2

#Number of the job in question
RUN=$4
export NJOB=$(($1+1+$RUN*10000))
echo "NJOB= "$NJOB

export energy=$3

#exporting eviroment variables necessary for the run
export MacroDir="/sampa/fcanedo/jewel-2.2.0"
export LD_LIBRARY_PATH="/sampa/fcanedo/lhapdf/lib"
export LHAPATH="/sampa/fcanedo/lhapdf/share/lhapdf/PDFsets"
export RIVET_ANALYSIS_PATH="/sampa/fcanedo/rivetanalises"
export LOG_FOLDER="/sampa/fcanedo/jewel-2.2.0/logs/pp_"$energy"_GeV/"
export LOG_FILE=$LOG_FOLDER"pp_"$tc"_$mdsfactor_."$NJOB".log"
export HEPMC_FOLDER="/sampa/fcanedo/jewel-2.2.0/hepmc/pp_"$energy"_GeV/"
export HEPMC_FILE=$HEPMC_FOLDER"pp_"$energy"_."$NJOB".hepmc"
export YODA_FOLDER="/sampa/fcanedo/results/histograms/pp_"$energy"_GeV/"
export YODA_GENERIC_FILE=$YODA_FOLDER"pp_"$energy"_."$NJOB".yoda"
export YODA_RAA_FILE=$YODA_FOLDER"raa.pp_"$energy"_."$NJOB".yoda"
export YODA_SHAPE_FILE=$YODA_FOLDER"shape.pp_"$energy"_."$NJOB".yoda"
export YODA_VN_FILE=$YODA_FOLDER"vn.pp_"$energy"_."$NJOB".yoda"
export YODA_ZHADRON_FILE=$YODA_FOLDER"zhadron.pp_"$energy"_."$NJOB".yoda"

if [ ! -e $YODA_FOLDER ] ; then
      mkdir -p $YODA_FOLDER
fi

if [ ! -e $LOG_FOLDER ] ; then
      mkdir -p $LOG_FOLDER
fi

if [ ! -e $HEPMC_FOLDER ] ; then
      mkdir -p $HEPMC_FOLDER
fi


chmod 777 $MacroDir/jewel-2.2.0-simple


source /cvmfs/alice.cern.ch/etc/login.sh

eval `alienv printenv VO_ALICE@GSL::v1.16-25`                                   
#eval `alienv printenv VO_ALICE@boost::v1.59.0-21`                               
#eval `alienv printenv VO_ALICE@HepMC::v2.06.09-13`                              
#eval `alienv printenv VO_ALICE@fastjet::v3.2.1_1.024-alice3-1`                  
#eval `alienv printenv VO_ALICE@YODA::v1.7.0-1`                                  
#eval `alienv printenv VO_ALICE@Rivet::2.6.0-alice1-1`                           
eval `alienv printenv VO_ALICE@Rivet::2.7.2-alice2-1`                           


#Decide whether to perform simulation or just analysis
if true ; then

      if [ -e $HEPMC_FILE".xz" ]; then
            rm $HEPMC_FILE".xz"
      fi

      cd $MacroDir/params/pp_"$energy"_GeV/
      cp pp.dat pp."$NJOB".dat

      sed -i "s/\(NJOB\s\).*/\1"$NJOB"/g" pp."$NJOB".dat
      sed -i "s/\(NEVENT\s\).*/\1"$NEVENT"/g" pp."$NJOB".dat
      sed -i "s|\(LOGFILE\s\).*|\1$LOG_FILE|g" pp."$NJOB".dat
      sed -i "s|\(HEPMCFILE\s\).*|\1$HEPMC_FILE|g" pp."$NJOB".dat
      sed -i "s/\(MEDIUMPARAMS\s.*pp\.\).*/\1"$NJOB".dat/g" pp."$NJOB".dat
      sed -i "s/\(XSECFILE\s.*pp_\).*/\1"$NJOB".dat/g" pp."$NJOB".dat
      sed -i "s/\(SPLITINTFILE\s.*pp_\).*/\1"$NJOB".dat/g" pp."$NJOB".dat
      sed -i "s/\(PDFFILE\s.*pp_\).*/\1"$NJOB".dat/g" pp."$NJOB".dat

      #cat "pp."$NJOB".dat"

      cd $MacroDir/medparams/pp_"$energy"_GeV/
      cp pp.dat pp."$NJOB".dat
      #cat pp.dat

      #cd $MacroDir/integrals
      #cp xsecs.dat xsecs_pp_$NJOB".dat"
      #cp splitint.dat splitint_pp_$NJOB".dat"
      #cp pdfs.dat pdfs_pp_$NJOB".dat"

      cd $MacroDir
      $MacroDir/jewel-2.2.0-vac $MacroDir/params/pp_2760_GeV/pp."$NJOB".dat

      cd $MacroDir/integrals
      rm xsecs_pp_$NJOB".dat"
      rm splitint_pp_$NJOB".dat"
      rm pdfs_pp_$NJOB".dat"

      cd $MacroDir/params/pp_"$energy"_GeV/
      rm pp."$NJOB".dat

      cd $MacroDir/medparams/pp_"$energy"_GeV/
      rm pp."$NJOB".dat

fi

if [ -e $HEPMC_FILE".xz" ]; then
      rm $HEPMC_FILE".xz"
fi

echo $HEPMC_FILE

export PSI=0

#rivet -a MC_GENERIC --ignore-beams -H $YODA_GENERIC_FILE $HEPMC_FILE
rivet  -a JET_SHAPE --ignore-beams -H $YODA_SHAPE_FILE $HEPMC_FILE
rivet  -a JET_RAA --ignore-beams -H $YODA_RAA_FILE $HEPMC_FILE
#rivet  -a JET_VN --ignore-beams -H $YODA_VN_FILE $HEPMC_FILE
#rivet -a Z_HADRON --ignore-beams -H $YODA_ZHADRON_FILE $HEPMC_FILE

#xz $HEPMC_FILE
