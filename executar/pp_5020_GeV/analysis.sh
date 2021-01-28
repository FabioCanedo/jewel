#!/bin/bash
cd /sampa/fcanedo/jewel-2.2.0/executar/pp_5020_GeV

BEGIN=$((100*$1))
END=$(($BEGIN+99))

for I in `seq $BEGIN $END` ; do

      bash pp.sh $I

done
