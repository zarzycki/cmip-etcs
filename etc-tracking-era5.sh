#!/bin/bash -l

############ USER OPTIONS #####################

## If using unstructured data, need connect file, otherwise empty
CONNECTFLAG="" 
DO_TRACKS=false
DO_EXTRACT=true
UQSTR=ERA5

############ MACHINE SPECIFIC AUTO-CONFIG #####################

  TEMPESTEXTREMESDIR=/glade/u/home/zarzycki/work/tempestextremes/
  PATHTOFILES=/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/
  PATHTOFILES2=/glade/u/home/zarzycki/rda/ds633.0/e5.oper.fc.sfc.meanflux/
  PATHTOFILES3=/glade/scratch/zarzycki/CMIPTMP/${UQSTR}/
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
  find ${PATHTOFILES} -name "e5.oper.an.sfc.128_151_msl.*nc" | grep .1986 | grep -v .2020 | sort -n >> $FILELISTNAME
  # Add any static file(s) to each line
  #$THISSED -e 's?$?;'"${TOPOFILE}"'?' -i $FILELISTNAME
  #cat $FILELISTNAME

  DCU_PSLNAME="MSL"    # Name of PSL on netcdf files
  DCU_PSLFOMAG=200.0   # Pressure "anomaly" required for PSL min to be kept
  DCU_PSLFODIST=6.0    # Pressure "anomaly" contour must lie within this distance
  DCU_MERGEDIST=6.0    # Merge two competing PSL minima within this radius
  SN_TRAERA5NGE=6.0     # Maximum distance (GC degrees) ETC can move in successive steps
  SN_TRAJMINLENGTH=10  # Min length of trajectory (nsteps)
  SN_TRAJMAXGAP=3      # Max gap within a trajectory (nsteps)

  STRDETECT="--verbosity 0 --timestride 6 --latname latitude --lonname longitude --out cyclones_tempest.${DATESTRING} --closedcontourcmd ${DCU_PSLNAME},${DCU_PSLFOMAG},${DCU_PSLFODIST},0 --mergedist ${DCU_MERGEDIST} --searchbymin ${DCU_PSLNAME} --outputcmd ${DCU_PSLNAME},min,0"
  echo $STRDETECT
  touch cyclones.${DATESTRING}
  mpiexec_mpt ${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STRDETECT} </dev/null

  cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
  rm cyclones_tempest.${DATESTRING}*

  # Stitch candidate cyclones together
  STRSTITCH="--range ${SN_TRAERA5NGE} --mintime 18h --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist 12.0 --threshold lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"
  ${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_fmt "lon,lat,slp" ${STRSTITCH}

  # Clean up leftover files
  rm cyclones.${DATESTRING}   #Delete candidate cyclone file
  rm ${FILELISTNAME}
  rm log*.txt

  ${TEMPESTEXTREMESDIR}/bin/HistogramNodes --in ${TRAJFILENAME} --nlat 36 --nlon 72 --iloncol 3 --ilatcol 4 --out "${PATHTOFILES3}/density_${UQSTR}.nc"
fi

############ FIELD EXTRACTION #####################
######## OK, do some node work now that we've tracked

if [ "$DO_EXTRACT" = true ] ; then

  touch $NODELISTNAME

  find ${PATHTOFILES2} -name "e5.oper.fc.sfc.meanflux.235_055_mtpr.ll025sc.*.nc" | grep .1986 | grep -v .2020 | sort -n >> $NODELISTNAME

  ## Create files and folder to hold filtered data
  cp $NODELISTNAME $NODELISTNAMEFILT
  $THISSED -i 's/mtpr_/mtpr_FILT_/g' $NODELISTNAMEFILT
  $THISSED -i 's?'${PATHTOFILES2//\?/\.}'?'${PATHTOFILES3}'?g' $NODELISTNAMEFILT
  mkdir -p ${PATHTOFILES3}

  ${TEMPESTEXTREMESDIR}/bin/NodeFileFilter --latname "latitude" --lonname "longitude" --in_nodefile ${TRAJFILENAME} --in_fmt "lon,lat,slp" --in_data_list ${NODELISTNAME} --out_data_list ${NODELISTNAMEFILT} --var "MTPR" --nearbyblobs "MTPR,2.0,>=,5.0e-5,50.0" --maskvar "mask" #--bydist 25.0    pr is kg/m2/s so pr/1000 => m/s


  pattern=${PATHTOFILES3}"/mtpr_FILT*${UQSTR}*.nc"
  dumfiles=( $pattern )
  ${TEMPESTEXTREMESDIR}/bin/NodeFileEditor --in_nodefile ${TRAJFILENAME} --in_data ${dumfiles[0]} --in_fmt "lon,lat,slp" --out_nodefile ${TRAJFILENAME}_strong.txt --out_fmt "lon,lat,slp" --col_filter "slp,<=,97000."

# in hypothetical world, where flatten time would go
#  for f in /glade/u/home/zarzycki/scratch/CMIPTMP/NESM3/pr_*nc
#  do
#    echo $f
#    ncl flatten-time.ncl 'filename="'${f}'"'
#  done

  ${TEMPESTEXTREMESDIR}/bin/NodeFileCompose --in_nodefile ${TRAJFILENAME}_strong.txt --in_fmt "lon,lat,slp" --in_data_list "${NODELISTNAMEFILT}" --var "MTPR" --max_time_delta "1h" --out_data "${PATHTOFILES3}/composite_${UQSTR}.nc" --dx 1.0 --resx 80 --op "mean" --snapshots
  
  # Clean up
  rm ${NODELISTNAME}
  rm ${NODELISTNAMEFILT}

fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime},${TRAJFILENAME}\n" >> timing.txt


