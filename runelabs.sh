#!/bin/bash
workdir="C:/Users/YK/Desktop/biot-res/"
# workdir="./biot-res/"  # use as a detault 
read -p "Type the name of simulation, followed by [enter]:" workname 
today=$(date "+%Y%m%d")
mkdir -p "$workdir$workname-$today"
make 
cp labs eParameters_IN.txt Parameters_IN.txt SedENV.IN Organisms.IN "$workdir$workname-$today"
$workdir$workname-$today/labs  "$workdir$workname-$today"