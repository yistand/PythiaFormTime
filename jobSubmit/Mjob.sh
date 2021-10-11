#!/bin/bash


#Array declaration
declare -a ptarr=(0 2 4 6 8 10 15 20 25 30 35 40 50 60 100)	#  first 0 for softqcd with pTHat_Max cut off

Msub=100

while [ $Msub -gt 0 ]
do
if [[ `condor_q | grep "Total for yili: 0 jobs"` ]]
then
for ((i=1; i<${#ptarr[@]}; i++))
do
echo submit -- $Msub -- pt$i
echo condor_submit ptmin="${ptarr[$(($i-1))]}" ptmax="${ptarr[$i]}"  job.jdl
condor_submit ptmin="${ptarr[$(($i-1))]}" ptmax="${ptarr[$i]}" job.jdl
done
let Msub-=1
fi

sleep 10


done
