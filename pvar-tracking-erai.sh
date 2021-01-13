#!/bin/bash -l

############ USER OPTIONS #####################

## If using unstructured data, need connect file, otherwise empty
CONNECTFLAG="" 
DO_TRACKS=true
UQSTR=ERAI-PVAR

############ MACHINE SPECIFIC AUTO-CONFIG #####################

TEMPESTEXTREMESDIR=/glade/u/home/zarzycki/work/tempestextremes/
PATHTOFILES=/glade/scratch/zarzycki/erain/
THISSED="sed"

############ TRACKER MECHANICS #####################
starttime=$(date -u +"%s")

DATESTRING=`date +"%s%N"`
OUTTRAJDIR="./trajs/"
TRAJFILENAME=${OUTTRAJDIR}/trajectories.txt.${UQSTR}
FILELISTNAME=filelist.txt.${DATESTRING}
NODELISTNAME=nodelist.txt.${DATESTRING}
NODELISTNAMEFILT=nodelistfilt.txt.${DATESTRING}

############ TRACK ETCs #####################

## Build trajectory files; the most time intensive part of this.
if [ "$DO_TRACKS" = true ] ; then

  mkdir -p $OUTTRAJDIR

  touch $FILELISTNAME

  #### Files
  find ${PATHTOFILES} -name "erain_*_pvar_one_side.nc" | sort -n >> $FILELISTNAME
  # Add any static file(s) to each line
  #$THISSED -e 's?$?;'"${TOPOFILE}"'?' -i $FILELISTNAME
  #cat $FILELISTNAME

  DCU_PSLNAME="pvar_dbl"    # Name of PSL on netcdf files
  DCU_PSLFOMAG=-200.0   # Pressure "anomaly" required for PSL min to be kept
  DCU_PSLFODIST=8.0    # Pressure "anomaly" contour must lie within this distance
  DCU_MERGEDIST=8.0    # Merge two competing PSL minima within this radius
  
  SN_TRAJRANGE=6.0     # Maximum distance (GC degrees) ETC can move in successive steps
  SN_TRAJMINLENGTH=10  # Min length of trajectory (nsteps)
  SN_TRAJMAXGAP=3      # Max gap within a trajectory (nsteps)

#--closedcontourcmd ${DCU_PSLNAME},${DCU_PSLFOMAG},${DCU_PSLFODIST},0   --minlat 25. --maxlat 70. 
  STRDETECT="--verbosity 0 --timestride 1 --latname lat --lonname lon --out cyclones_tempest.${DATESTRING} --closedcontourcmd ${DCU_PSLNAME},${DCU_PSLFOMAG},${DCU_PSLFODIST},0 --mergedist ${DCU_MERGEDIST} --searchbythreshold >0 --searchbymax ${DCU_PSLNAME} --outputcmd ${DCU_PSLNAME},max,0"
  echo $STRDETECT
  touch cyclones.${DATESTRING}
  mpiexec_mpt ${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STRDETECT} </dev/null

  cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
  rm cyclones_tempest.${DATESTRING}*

  # Stitch candidate cyclones together
  #--threshold lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1
  STRSTITCH="--range ${SN_TRAJRANGE} --mintime 18h --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist 12.0"
  ${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_fmt "lon,lat,slp" ${STRSTITCH}

  # Clean up leftover files
  rm cyclones.${DATESTRING}   #Delete candidate cyclone file
  rm ${FILELISTNAME}
  rm log*.txt

  ${TEMPESTEXTREMESDIR}/bin/HistogramNodes --in ${TRAJFILENAME} --nlat 36 --nlon 72 --iloncol 3 --ilatcol 4 --out "${PATHTOFILES}/density_${UQSTR}.nc"
fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime},${TRAJFILENAME}\n" >> timing.txt


