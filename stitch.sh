#!/bin/bash

CONNECTFLAG="" 
DO_TRACKS=true
DO_EXTRACT=false
ENSEMBLEMEMBER=r1i1p1f1

UQSTR="ERAI"
PARENTSTR=""
GRIDDATE=""
STYR="-1"

DATESTRING="1612378751510351441"
TRAJFILENAME="traj.out"

############ MACHINE SPECIFIC AUTO-CONFIG #####################

TEMPESTEXTREMESDIR=/glade/u/home/zarzycki/work/tempestextremes/
PATHTOFILES=/glade/scratch/zarzycki/h1files/${UQSTR}/
PATHTOFILES3=/glade/scratch/zarzycki/CMIPTMP/${UQSTR}/
TOPOFILE=/glade/work/zarzycki/topo-files/CMIP6/${UQSTR}.topo.nc
TOPOORIGDIR=/glade/work/zarzycki/cam_tools/hires-topo/
THISSED="sed"

############ MACHINE SPECIFIC AUTO-CONFIG #####################

SN_TRAJRANGE=6.0     # Maximum distance (GC degrees) ETC can move in successive steps
SN_TRAJMINLENGTH=12  # Min length of trajectory (nsteps)
SN_TRAJMAXGAP=0      # Max gap within a trajectory (nsteps)
SN_MINENDPOINT=12.0  # Min travel distance

# Stitch candidate cyclones together
STRSTITCH="--range ${SN_TRAJRANGE} --mintime 60h --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist ${SN_MINENDPOINT} --threshold lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"
${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_fmt "lon,lat,slp,phis" ${STRSTITCH}