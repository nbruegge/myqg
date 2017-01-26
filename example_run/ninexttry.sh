#!/bin/bash
dir="../run_$1"
if [ -d ${dir} ]; then
    echo Warning the directory ${dir} does already exist!
    exit
fi
mkdir $dir
cp ni* setup_*.f90 $dir
