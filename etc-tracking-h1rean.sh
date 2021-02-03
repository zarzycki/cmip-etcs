#!/bin/bash -l

#PBS -N tempest.par
#PBS -A UNSB0017 
#PBS -l walltime=0:59:00
#PBS -q premium
#PBS -j oe
#PBS -l select=4:ncpus=36:mpiprocs=36

############ USER OPTIONS #####################

## If using unstructured data, need connect file, otherwise empty
CONNECTFLAG="" 
DO_TRACKS=true
DO_EXTRACT=false
ENSEMBLEMEMBER=r1i1p1f1

UQSTR="NARR"
PARENTSTR=""
GRIDDATE=""
STYR="-1"

############ MACHINE SPECIFIC AUTO-CONFIG #####################

if [[ $HOSTNAME = cheyenne* ]]; then
  TEMPESTEXTREMESDIR=/glade/u/home/zarzycki/work/tempestextremes/
  PATHTOFILES=/glade/scratch/zarzycki/h1files/${UQSTR}/
  PATHTOFILES3=/glade/scratch/zarzycki/CMIPTMP/${UQSTR}/
  TOPOFILE=/glade/work/zarzycki/topo-files/CMIP6/${UQSTR}.topo.nc
  TOPOORIGDIR=/glade/work/zarzycki/cam_tools/hires-topo/
  THISSED="sed"
  DNCOMMAND="mpiexec_mpt"
else
  echo "Can't figure out hostname, exiting"
  exit
fi

############ TRACKER MECHANICS #####################
starttime=$(date -u +"%s")

DATESTRING=`date +"%s%N"`
OUTTRAJDIR="./trajs/"
TRAJFILENAME=${OUTTRAJDIR}/trajectories.txt.${ENSEMBLEMEMBER}.${UQSTR}
FILELISTNAME=filelist.txt.${DATESTRING}
NODELISTNAME=nodelist.txt.${DATESTRING}
NODELISTNAMEFILT=nodelistfilt.txt.${DATESTRING}

########### SET UP FOLDERS #####################

mkdir -p ${PATHTOFILES3}

############ TRACK ETCs #####################

## Build trajectory files; the most time intensive part of this.
if [ "$DO_TRACKS" = true ] ; then

  mkdir -p $OUTTRAJDIR

  touch $FILELISTNAME

  #### Files
  find ${PATHTOFILES} -name "*.h1.*.nc" | sort -n >> $FILELISTNAME
  # Add any static file(s) to each line
  $THISSED -e 's?$?;'"${TOPOFILE}"'?' -i $FILELISTNAME
  cat $FILELISTNAME

  ## Check if topo file exists, if not generate one using the first PSL file
  if [ ! -f ${TOPOFILE} ]; then
    echo "Topo file ${TOPOFILE} not found!"
    FIRSTPSLFILE=$(echo $(head -n 1 ${FILELISTNAME}) | awk -F';' '{print $1}')
    export TOPOORIGDIR="${TOPOORIGDIR}"  # Send topo dir to shell var for use inside NCL
    ncl data-process/make_topo_file.ncl 'trackerFileName="'$FIRSTPSLFILE'"' 'outFileName="'${TOPOFILE}'"'
  else
    echo "Topo file ${TOPOFILE} already exists."
  fi

  DCU_PSLNAME="PSL"    # Name of PSL on netcdf files
  DCU_PSLFOMAG=200.0   # Pressure "anomaly" required for PSL min to be kept
  DCU_PSLFODIST=6.0    # Pressure "anomaly" contour must lie within this distance
  DCU_MERGEDIST=6.0    # Merge two competing PSL minima within this radius
  SN_TRAJRANGE=6.0     # Maximum distance (GC degrees) ETC can move in successive steps
  SN_TRAJMINLENGTH=10  # Min length of trajectory (nsteps)
  SN_TRAJMAXGAP=3      # Max gap within a trajectory (nsteps)
  SN_MINENDPOINT=12.0  # Min travel distance

  STRDETECT="--verbosity 0 --timestride 1 ${CONNECTFLAG} --out cyclones_tempest.${DATESTRING} --closedcontourcmd ${DCU_PSLNAME},${DCU_PSLFOMAG},${DCU_PSLFODIST},0 --mergedist ${DCU_MERGEDIST} --searchbymin ${DCU_PSLNAME} --outputcmd ${DCU_PSLNAME},min,0;PHIS,max,0"
  echo $STRDETECT
  touch cyclones.${DATESTRING}
  ${DNCOMMAND} ${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STRDETECT} </dev/null
  cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
  rm cyclones_tempest.${DATESTRING}*

  ### Special filtering by time to account for fact that some models have longer PSL datasets than PRECT or PRECSN
  ## Create a temp cyclones file only including dates from STYR onwards, then overwrite full cyclones file before SN
  if [ $STYR \> 0 ]; then
    ${THISSED} -n '/^'${STYR}'/,$p' cyclones.${DATESTRING} > cycfilter.tmp.${DATESTRING}.txt;
    mv cycfilter.tmp.${DATESTRING}.txt cyclones.${DATESTRING}
  fi;

  # Stitch candidate cyclones together
  STRSTITCH="--range ${SN_TRAJRANGE} --mintime 60h --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist ${SN_MINENDPOINT} --threshold lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"
  ${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_fmt "lon,lat,slp,phis" ${STRSTITCH}

  # Clean up leftover files
  rm cyclones.${DATESTRING}   #Delete candidate cyclone file
  rm ${FILELISTNAME}
  rm log*.txt

  ${TEMPESTEXTREMESDIR}/bin/HistogramNodes --in ${TRAJFILENAME} --nlat 36 --nlon 72 --iloncol 3 --ilatcol 4 --out "${PATHTOFILES3}/density_${ENSEMBLEMEMBER}_${UQSTR}.nc"
fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime},${TRAJFILENAME}\n" >> timing.txt
