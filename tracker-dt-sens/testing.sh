#!/bin/bash -l

DCU_PSLNAME="MSL"    # Name of PSL on netcdf files
DCU_PSLFOMAG=200.0   # Pressure "anomaly" required for PSL min to be kept
DCU_PSLFODIST=6.0    # Pressure "anomaly" contour must lie within this distance
DCU_MERGEDIST=9.01   # Merge two competing PSL minima within this radius

FREQ="01"
SN_TRAJRANGE=$1
SN_TRAJMINLENGTH=$2
SN_TRAJMAXGAP=$3
SN_MINENDPOINT=12.0
SN_MINTIME=$4
DATESTRING=1708552749915423205

#FREQ="06"
#SN_TRAJRANGE=9.01   # Maximum distance (GC degrees) ETC can move in successive steps
#SN_TRAJMINLENGTH=4 # Min length of trajectory (nsteps)
#SN_TRAJMAXGAP=1      # Max gap within a trajectory (nsteps)
#SN_MINENDPOINT=12.0  # Min travel distance
#SN_MINTIME="24h"
#DATESTRING=1708624992021034724 # 6hourly

UQSTR=$FREQ-RANGE_$SN_TRAJRANGE-MINLEN_$SN_TRAJMINLENGTH-MAXGAP_$SN_TRAJMAXGAP-MINEND_$SN_MINENDPOINT-MINTIME_$SN_MINTIME

TEMPESTEXTREMESDIR=/glade/u/home/zarzycki/work/tempestextremes_noMPI/
TRAJFILENAME="traj-${UQSTR}.txt"

STRSTITCH="--range ${SN_TRAJRANGE} --mintime ${SN_MINTIME} --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist ${SN_MINENDPOINT}"
${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_fmt "lon,lat,psl,phis" ${STRSTITCH}

${TEMPESTEXTREMESDIR}/bin/HistogramNodes --in ${TRAJFILENAME} --nlat 90 --nlon 180 --iloncol 3 --ilatcol 4 --out "histogram-${UQSTR}.nc"

