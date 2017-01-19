#!/bin/bash

# define setup file (relative to path_model)
setup_file="setup_double_gyre.f90"
path_model="/Users/nbruggemann/syn_dropbox/sync/Promotion/src/myqg/"
path_act=`pwd`

## get filename without path and ending
##setup=`basename ${setup_file} .f90`

# cleaning up old run
rm ferret.jnl* *.data *.meta pybintocdf.cdf

# compile the model and copy executable to this path
cd $path_model
rm ./tmp_setup.f90
cp -r ${path_act}/${setup_file} ./tmp_setup.f90
make setup="tmp_setup"
cd $path_act
cp $path_model/myqg.out ./

# run the model
echo "start model run at `date +%Y-%m-%d_%H:%M:%S`" >> timing.txt
./myqg.out
echo "end   model run at `date +%Y-%m-%d_%H:%M:%S`" >> timing.txt

# make netcdf file
#./nipybintocdf.py
