#!/bin/bash -l

#PBS -l nodes=1:ppn=16
#PBS -l walltime=8:00:00
#PBS -A cmz5202_a_g_sc_default

#ncks -d lat,,,8 -d lon,,,8 /storage/home/cmz5202/work/cam_tools/hires-topo/2deg_cesm_topo_latlon.nc tmp.nc
#ncremap -a patch -i tmp.nc -o WRF2.nc -d ~/group/TGW/pgw/wrfout_d01_1995-06-04_00\:00\:00_3hourly.aux.nc

cd /storage/home/cmz5202/work/sw/cmip-etcs/data-process

TOPOFILE=WRF.nc
FIRSTPSLFILE=/storage/home/cmz5202/group/TGW/pgw/wrfout_d01_1995-06-04_00:00:00_3hourly.aux.nc
TOPOORIGDIR=/storage/home/cmz5202/work/cam_tools/hires-topo/

export TOPOORIGDIR="${TOPOORIGDIR}"
ncl make_topo_file.ncl 'trackerFileName="'$FIRSTPSLFILE'"' 'outFileName="'${TOPOFILE}'"'

