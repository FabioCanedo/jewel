#!/bin/bash

#Usage:
#$1=NJOB
#$2=NEVENT
#$3=energy
#$4=mdscalefactor
#$5=minimum centrality
#$6=maximum centrality
#$7=critical temperature
#$8=RUN


NEVENT=$2

#Number of the job in question
RUN=$8
export NJOB=$(($1+1+$RUN*10000))
echo "NJOB= "$NJOB

export energy=$3
export mdsfactort=$4
export mincent=$5
export maxcent=$6
export tc=$7

#exporting eviroment variables necessary for the run
export FileDir="/sampa/fcanedo/IC/mymedf/PbPb_"$energy"_GeV/"$mincent"_"$maxcent
export MacroDir="/sampa/fcanedo/jewel-2.2.0"
export LD_LIBRARY_PATH="/sampa/fcanedo/lhapdf/lib"
export LHAPATH="/sampa/fcanedo/lhapdf/share/lhapdf/PDFsets"
export RIVET_ANALYSIS_PATH="/sampa/fcanedo/rivetanalises"
export LOG_FOLDER="/sampa/fcanedo/jewel-2.2.0/logs/PbPb_"$energy"_GeV/cent"$mincent"_"$maxcent"/"
export LOG_FILE=$LOG_FOLDER"read_"$tc"_$mdsfactor_."$NJOB".log"
export HEPMC_FOLDER="/sampa/fcanedo/jewel-2.2.0/hepmc/PbPb_"$energy"_GeV/cent"$mincent"_"$maxcent"/"
export HEPMC_FILE=$HEPMC_FOLDER"read_"$tc"_$energy_."$NJOB".hepmc"
export YODA_FOLDER="/sampa/fcanedo/results/histograms/PbPb_"$energy"_GeV/cent"$mincent"_"$maxcent"/"
export YODA_GENERIC_FILE=$YODA_FOLDER"read_"$tc"_"$energy"_."$NJOB".yoda"
export YODA_RAA_FILE=$YODA_FOLDER"raa.read_"$tc"_"$energy"_."$NJOB".yoda"
export YODA_SHAPE_FILE=$YODA_FOLDER"shape.read_"$tc"_"$energy"_."$NJOB".yoda"
export YODA_VN_FILE=$YODA_FOLDER"vn.read_"$tc"_"$energy"_."$NJOB".yoda"
export YODA_ZHADRON_FILE=$YODA_FOLDER"zhadron.read_"$tc"_"$energy"_."$NJOB".yoda"

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

export MED_FILE_NUMBER=$(($1+1))

#Decide whether to perform simulation or just analysis
if true ; then

      echo "$FileDir"/"$MED_FILE_NUMBER".dat.xz
      if [ -e "$FileDir"/"$MED_FILE_NUMBER".dat.xz ]; then
            xz -d "$FileDir"/"$MED_FILE_NUMBER".dat.xz
            echo "Medfile exists"
      else
            echo "Medfile does not exist"
      fi

      if [ -e $HEPMC_FILE".xz" ]; then
            rm $HEPMC_FILE".xz"
      fi

      cd $MacroDir/params/PbPb_"$energy"_GeV/
      cp read.dat read."$NJOB".dat

      sed -i "s/\(NJOB\s\).*/\1"$NJOB"/g" read."$NJOB".dat
      sed -i "s/\(NEVENT\s\).*/\1"$NEVENT"/g" read."$NJOB".dat
      sed -i "s|\(LOGFILE\s\).*|\1$LOG_FILE|g" read."$NJOB".dat
      sed -i "s|\(HEPMCFILE\s\).*|\1$HEPMC_FILE|g" read."$NJOB".dat
      sed -i "s/\(MEDIUMPARAMS\s.*read\.\).*/\1"$NJOB".dat/g" read."$NJOB".dat
      sed -i "s/\(XSECFILE\s.*read_\).*/\1"$NJOB".dat/g" read."$NJOB".dat
      sed -i "s/\(SPLITINTFILE\s.*read_\).*/\1"$NJOB".dat/g" read."$NJOB".dat
      sed -i "s/\(PDFFILE\s.*read_\).*/\1"$NJOB".dat/g" read."$NJOB".dat

      cd $MacroDir/medparams/PbPb_"$energy"_GeV/
      cp read.dat read."$NJOB".dat
      #cat read.dat

      sed -i "s/\(MEDFILE\s.*\/\)[0-9]*\(\.dat\)/\1"$MED_FILE_NUMBER"\2/g" read."$NJOB".dat
      sed -i "s/\(MDSCALEFAC\s\).*/\1$mdsfactord0/g" read."$NJOB".dat

      sed -i "s|\(CENTRMIN\s\).*|\1"$mincent".|g" read."$NJOB".dat
      sed -i "s|\(CENTRMAX\s\).*|\1"$maxcent".|g" read."$NJOB".dat
      sed -i "s|\(TC\s\).*|\1"$tc"d0|g" read."$NJOB".dat

      #cd $MacroDir/integrals
      #cp xsecs.dat xsecs_read_$NJOB".dat"
      #cp splitint.dat splitint_read_$NJOB".dat"
      #cp pdfs.dat pdfs_read_$NJOB".dat"

      cd $MacroDir
      $MacroDir/jewel-2.2.0-reader $MacroDir/params/PbPb_2760_GeV/read."$NJOB".dat

      cd $MacroDir/integrals
      rm xsecs_read_$NJOB".dat"
      rm splitint_read_$NJOB".dat"
      rm pdfs_read_$NJOB".dat"

      cd $MacroDir/params/PbPb_"$energy"_GeV/
      rm read."$NJOB".dat

      cd $MacroDir/medparams/PbPb_"$energy"_GeV/
      rm read."$NJOB".dat

      if [ ! -e "$FileDir"/"$MED_FILE_NUMBER".dat.xz ]; then
            xz "$FileDir"/"$MED_FILE_NUMBER".dat
      fi

fi

export PSI=0.0
export PSI2=0.0
export PSI3=0.0
export PSI4=0.0

if [ -e $HEPMC_FILE".xz" ]; then
      rm $HEPMC_FILE".xz"
fi

echo $HEPMC_FILE

#rivet -a MC_GENERIC --ignore-beams -H $YODA_GENERIC_FILE $HEPMC_FILE
rivet  -a JET_SHAPE --ignore-beams -H $YODA_SHAPE_FILE $HEPMC_FILE
rivet  -a JET_RAA --ignore-beams -H $YODA_RAA_FILE $HEPMC_FILE
rivet  -a JET_VN --ignore-beams -H $YODA_VN_FILE $HEPMC_FILE
#rivet -a Z_HADRON --ignore-beams -H $YODA_ZHADRON_FILE $HEPMC_FILE

xz $HEPMC_FILE
