#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -N qg_wind_gyre
#PBS -o log-$PBS_JOBID.out
#PBS -e log-$PBS_JOBID.err

if [ "x$PBS_O_WORKDIR" != "x" ]
then
  cd $PBS_O_WORKDIR
fi

# define setup file (relative to path_model)
setup_file="setup_double_gyre.f90"
#path_model="/home/nbruggemann/src/myqg/"
path_model="/Users/nbruggemann/syn_dropbox/sync/Promotion/src/myqg/"
path_act=`pwd`

## get filename without path and ending
##setup=`basename ${setup_file} .f90`

# replace path_data variable in setup_file
sed -i -e s,"  path_data.*","  path_data = \"${path_act}/\"",g ${setup_file}

# cleaning up old run
#rm ferret.jnl* *.data *.meta pybintocdf.cdf
./niclear.sh

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
./nipybintocdf.py
