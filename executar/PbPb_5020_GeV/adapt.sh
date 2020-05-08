#!/bin/bash

for I in `seq 1 10`
do

FILENUMBER=`cat ~/IC/vusphydro/PbPb_5000_GeV/first_filelist.txt | head -n $I | tail -n 1`

sed -i "2s/\(.*\)/export FileDir=\/sampa\/fcanedo\/IC\/vusphydro\/PbPb_5000_GeV\n\1/g" vusp.$I.sh

sed -i "7s/jewel.*jakiEventFiles\/[0-9]*\.dat/IC\/vusphydro\/PbPb_5000_GeV\/$FILENUMBER"".dat.xz/g" vusp.$I.sh
sed -i "8s/jewel.*jakiEventFiles\/[0-9]*\.dat\.xz/IC\/vusphydro\/PbPb_5000_GeV\/$FILENUMBER"".dat.xz/g" vusp.$I.sh
sed -i "12s/jewel.*jakiEventFiles\/[0-9]*\.dat\.xz/IC\/vusphydro\/PbPb_5000_GeV\/$FILENUMBER"".dat.xz/g" vusp.$I.sh
sed -i "13s/jewel.*jakiEventFiles\/[0-9]*\.dat/IC\/vusphydro\/PbPb_5000_GeV\/$FILENUMBER"".dat.xz/g" vusp.$I.sh

sed -i "11s/\(params\)/\1\/PbPb_5000_GeV/g" vusp.$I.sh

sed -i "29s/\(.*\)/#\1\nrivet -a Z_HADRON --ignore-beams -H \/sampa\/fcanedo\/results\/yoda\/PbPb_5000_GeV\/cent0_10\/z_hadron\/vusp.$I.yoda \/sampa\/fcanedo\/yoda\/PbPb_5000_GeV\/cent0_10\/vusp.$I.hepmc/g" vusp.$I.sh

sed -i "7s/jewel.*jakiEventFiles\/[0-9]*\.dat/IC\/vusphydro\/PbPb_5000_GeV\/$FILENUMBER"".dat/g" ../../medparams/PbPb_5000_GeV/vusp.$I.dat

sed -i "6s/\(.*\)/\1\nPROCESS PPZJ/g" ../../params/PbPb_5000_GeV/vusp.$I.dat
sed -i "6s/\(.*\)/\1\nSQRTS 5020/g" ../../params/PbPb_5000_GeV/vusp.$I.dat
sed -i "3s/\(LOGFILE.*logs\)/\1\/PbPb_5020_GeV/g" ../../params/PbPb_5000_GeV/vusp.$I.dat
sed -i "4s/\(HEPMCFILE.*hepmc\)\/PbPb_2760/\1\/PbPb_5020_GeV/g" ../../params/PbPb_5000_GeV/vusp.$I.dat
sed -i "s/\(MEDIUMPARAMS.*medparams\/\).*\(vusp\)/\1PbPb_5000_GeV\/\2/g" ../../params/PbPb_5000_GeV/vusp.$I.dat

done
