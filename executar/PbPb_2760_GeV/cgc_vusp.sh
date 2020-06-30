#!/bin/bash

#Number of the job in question
export NJOB=5
echo "NJOB= "$NJOB

#exporting enviroment variables necessary for the run
export FileDir="/sampa/fcanedo/IC/vusphydro/cgc/profiles/0-5/"
export MacroDir="/sampa/fcanedo/jewel-2.2.0"
export LD_LIBRARY_PATH="/sampa/fcanedo/lhapdf/lib"
export LHAPATH="/sampa/fcanedo/lhapdf/share/lhapdf/PDFsets"
export RIVET_ANALYSIS_PATH="/sampa/fcanedo/rivetanalises"
export HEPMC_FILE="/sampa/fcanedo/jewel-2.2.0/hepmc/PbPb_2760_GeV/cent0_10/cgc_vusp."$NJOB".hepmc"
export YODA_GENERIC_FILE="/sampa/fcanedo/results/histograms/PbPb_2760_GeV/cent0_10/cgc_vusp/cgc_vusp."$NJOB".yoda"
export YODA_SHAPE_FILE="/sampa/fcanedo/results/histograms/PbPb_2760_GeV/cent0_10/shape_cgc_vusp/shape.cgc_vusp."$NJOB".yoda"
export YODA_VN_FILE="/sampa/fcanedo/results/histograms/PbPb_2760_GeV/cent0_10/vn_cgc_vusp/vn.cgc_vusp."$NJOB".yoda"
export YODA_ZHADRON_FILE="/sampa/fcanedo/results/histograms/PbPb_2760_GeV/cent0_10/z_hadron/cgc_vusp."$NJOB".yoda"

chmod 777 $MacroDir/jewel-2.2.0-simple

#Loading libraries and software
source /cvmfs/alice.cern.ch/etc/login.sh

eval `alienv printenv VO_ALICE@GSL::v1.16-25`                                   
eval `alienv printenv VO_ALICE@boost::v1.59.0-21`                               
eval `alienv printenv VO_ALICE@HepMC::v2.06.09-13`                              
eval `alienv printenv VO_ALICE@fastjet::v3.2.1_1.024-alice3-1`                  
eval `alienv printenv VO_ALICE@YODA::v1.7.0-1`                                  
eval `alienv printenv VO_ALICE@Rivet::2.6.0-alice1-1`                           

#Number used for the name of the medium file, not necessary for my cgc version
#export MED_FILE_NUMBER=`cat $FileDir/filelist.txt | head -n $NJOB | tail -n 1`
#echo "MED FILE NUMBER= "$MED_FILE_NUMBER

#This block was necessary for the medium conversion from root version(Caio), I already did it
#cd $FileDir
#sed -i "13s/\(BEGIN=\).*/\1$NJOB/g" reader.sh
#sed -i "14s/\(END=\).*/\1$NJOB/g" reader.sh
#./reader.sh

#Unpacking the compressed medium file
if [ -e "$FileDir"/"$NJOB".dat.xz ]; then
xz -d "$FileDir"/"$NJOB".dat.xz
echo "Medfile exists"
else
echo "Medfile does not exist"
fi

#Creating copies of the parameter files and making approriate changes for the Job number
cd $MacroDir/params/PbPb_2760_GeV/
cp cgc_vusp.dat cgc_vusp."$NJOB".dat

sed -i "s/\(NJOB\s\).*/\1"$NJOB"/g" cgc_vusp."$NJOB".dat
sed -i "s/\(LOGFILE\s.*cgc_vusp\.\).*/\1"$NJOB".log/g" cgc_vusp."$NJOB".dat
sed -i "s/\(HEPMCFILE\s.*cgc_vusp\.\).*/\1"$NJOB".hepmc/g" cgc_vusp."$NJOB".dat
sed -i "s/\(MEDIUMPARAMS\s.*cgc_vusp\.\).*/\1"$NJOB".dat/g" cgc_vusp."$NJOB".dat
sed -i "s/\(XSECFILE \s.*cgc_vusp_\).*/\1"$NJOB".dat/g" cgc_vusp."$NJOB".dat
sed -i "s/\(SPLITINTFILE\s.*cgc_vusp_\).*/\1"$NJOB".dat/g" cgc_vusp."$NJOB".dat
sed -i "s/\(PDFFILE\s.*cgc_vusp_\).*/\1"$NJOB".dat/g" cgc_vusp."$NJOB".dat

cd $MacroDir/medparams/PbPb_2760_GeV/
cp cgc_vusp.dat cgc_vusp."$NJOB".dat

sed -i "s/\(MEDFILE\s.*\/\)[0-9]*\(\.dat\)/\1"$NJOB"\2/g" cgc_vusp."$NJOB".dat

cd $MacroDir/integrals
cp xsecs.dat xsecs_cgc_trento_$NJOB".dat"
cp splitint.dat splitint_cgc_trento_$NJOB".dat"
cp pdfs.dat pdfs_cgc_trento_$NJOB".dat"

#Executing JEWEL
cd $MacroDir
$MacroDir/jewel-2.2.0-reader $MacroDir/params/PbPb_2760_GeV/cgc_vusp."$NJOB".dat

cd $MacroDir/integrals
rm xsecs_cgc_vusp_$NJOB".dat"
rm splitint_cgc_vusp_$NJOB".dat"
rm pdfs_cgc_vusp_$NJOB".dat"

#Cleaning the specific parameter files
cd $MacroDir/params/PbPb_2760_GeV/
rm cgc_vusp."$NJOB".dat

cd $MacroDir/medparams/PbPb_2760_GeV/
rm cgc_vusp."$NJOB".dat

#Compressing the medium files back
if [ ! -e "$FileDir"/"$NJOB".dat.xz ]; then
echo ""
xz "$FileDir"/"$NJOB".dat
fi



#let line=($NJOB%100+1)
line=`cat /sampa/fcanedo/IC/vusphydro/cgc/filelist_0_5.txt | head -n $NJOB | tail -n 1`
export PSI2=`cat /sampa/fcanedo/IC/vusphydro/cgc/eta=0.11_all_pt0.2-3.dat | head -n $line | tail -n 1 | sed "s/\s\{1,\}/ /g" | cut -d ' ' -f 5`
export PSI3=`cat /sampa/fcanedo/IC/vusphydro/cgc/eta=0.11_all_pt0.2-3.dat | head -n $line | tail -n 1 | sed "s/\s\{1,\}/ /g" | cut -d ' ' -f 7`
export PSI4=`cat /sampa/fcanedo/IC/vusphydro/cgc/eta=0.11_all_pt0.2-3.dat | head -n $line | tail -n 1 | sed "s/\s\{1,\}/ /g" | cut -d ' ' -f 9`
export PSI=$PSI2
echo $PSI

#Executing the analysis of the events
#rivet -a MC_GENERIC --ignore-beams -H $YODA_GENERIC_FILE $HEPMC_FILE
echo $YODA_SHAPE_FILE
echo $HEPMC_FILE
rivet  -a JET_SHAPE --ignore-beams -H $YODA_SHAPE_FILE $HEPMC_FILE
rivet  -a JET_VN --ignore-beams -H $YODA_VN_FILE $HEPMC_FILE
#rivet -a Z_HADRON --ignore-beams -H $YODA_ZHADRON_FILE $HEPMC_FILE
