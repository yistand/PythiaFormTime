#!/bin/bash


echo "This job is running on $(hostname)."
date

source /hep/home/yili/.bash_profile

npythia=$1
min=$2
max=$3
id=$4
nlbt=$5

## The first way to caculate Formation time
#cd ../pythiaFT
#echo $PWD
#
#echo
#make FormTime
#echo ./FormTime $npythia $min $max $id 
#./FormTime $npythia $min $max $id 
#
#echo
#
#
#
#cd ../JETSCAPE
#echo $PWD
#
#echo
#echo ./run-patch $nlbt $min $max $id
#./run-patch $nlbt $min $max $id
#echo
#
#
#cd ../JETSCAPE-noFT
#echo $PWD
#
#echo
#echo ./run-patch $nlbt $min $max $id
#./run-patch $nlbt $min $max $id
#echo


## The second way to caculate Formation time
#cd ../pythiaFT
#echo $PWD
#
#echo
#make FormTime_SimpleIF
#echo ./FormTime_SimpleIF $npythia $min $max $id 
#./FormTime_SimpleIF $npythia $min $max $id 
#
#echo
#
#
cd ../JETSCAPE
echo $PWD

echo
echo ./run-patch-SimpleIF $nlbt $min $max $id
./run-patch-SimpleIF $nlbt $min $max $id
echo




echo "Job Done!"
date
