#!/bin/bash

rm composite_CESM1-LENS.nc
cp composite_U.nc composite_CESM1-LENS.nc
ncks -A composite_V.nc composite_CESM1-LENS.nc
ncks -A composite_Q.nc composite_CESM1-LENS.nc
ncks -A composite_T.nc composite_CESM1-LENS.nc
ncks -A composite_PS.nc composite_CESM1-LENS.nc
ncks -A composite_PRECT.nc composite_CESM1-LENS.nc
ncks -A composite_PSL.nc composite_CESM1-LENS.nc
ncks -A composite_TREFHT.nc composite_CESM1-LENS.nc
ncks -A composite_QREFHT.nc composite_CESM1-LENS.nc
ncks -A -v hyai,hybi,hyam,hybm,lev,ilev /glade/u/home/zarzycki/betacast/atm_to_cam/templates/L30template.nc composite_CESM1-LENS.nc 
ncrename -v U\(:\),U composite_CESM1-LENS.nc
ncrename -v V\(:\),V composite_CESM1-LENS.nc
ncrename -v T\(:\),T composite_CESM1-LENS.nc
ncrename -v Q\(:\),Q composite_CESM1-LENS.nc

