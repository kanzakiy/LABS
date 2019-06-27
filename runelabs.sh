#!/bin/bash
workdir="C:/Users/YK/Desktop/biot-res/"
read -p "Type the name of simulation, followed by [ENTER]:" workname 
today=$(date "+%Y%m%d")
# echo $today
# echo $workdir
# echo $workname
mkdir -p "$workdir$workname-$today"
# make clean
make 
cp labs.exe "$workdir$workname-$today"
cp eParameters_IN.txt "$workdir$workname-$today"
cp Parameters_IN.txt "$workdir$workname-$today"
cp SedENV.IN  "$workdir$workname-$today"
cp Organisms.IN  "$workdir$workname-$today"
cd "$workdir$workname-$today"
$workdir$workname-$today/labs.exe  "$workdir$workname-$today"