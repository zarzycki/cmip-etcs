#!/bin/bash -l

############ USER OPTIONS #####################

## If using unstructured data, need connect file, otherwise empty
CONNECTFLAG="" 
DO_TRACKS=true
DO_EXTRACT=false
ENSEMBLEMEMBER=r1i1p1f1   #r1, r2, etc.

## Unique string

#UQSTR="GFDL-CM4"
#PARENTSTR="NOAA-GFDL"
#GRIDDATE="gr1/v20180701"
#STYR="1960"

#UQSTR="GISS-E2-1-G"
#PARENTSTR="NASA-GISS"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="MPI-ESM1-2-LR"
#PARENTSTR="MPI-M"
#GRIDDATE="gn/v????????"
#STYR="-1"

## CMZ: issues with snow files
#UQSTR="IPSL-CM6A-LR"
#PARENTSTR="IPSL"
#GRIDDATE="gr/v????????"
#STYR="-1"

#UQSTR="ACCESS-ESM1-5"
#PARENTSTR="CSIRO"
#GRIDDATE="gn/v????????"
#STYR="1960"

#UQSTR="MIROC6"
#PARENTSTR="MIROC"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="MRI-ESM2-0"
#PARENTSTR="MRI"
#GRIDDATE="gn/v????????"
#STYR="1960"

## CMZ: issues with time coordinate in pr/prsn files
#UQSTR="BCC-CSM2-MR"
#PARENTSTR="BCC"
#GRIDDATE="gn/v????????"
#STYR="1950"

# CMZ: need to flatten time
#UQSTR="NESM3"
#PARENTSTR="NUIST"
#GRIDDATE="gn/v????????"
#STYR="1960"

#UQSTR="SAM0-UNICON"
#PARENTSTR="SNU"
#GRIDDATE="gn/v????????"
#STYR="1950"

#UQSTR="NorESM2-LM"
#PARENTSTR="NCC"
#GRIDDATE="gn/v????????"
#STYR="-1"

UQSTR="AWI-ESM-1-1-LR"
PARENTSTR="AWI"
GRIDDATE="gn/v????????"
STYR="-1"

############ MACHINE SPECIFIC AUTO-CONFIG #####################

if [[ $HOSTNAME = cheyenne* ]]; then
  TEMPESTEXTREMESDIR=/glade/u/home/zarzycki/work/tempestextremes_noMPI/
  PATHTOFILES=/glade/collections/cmip/CMIP6/CMIP/${PARENTSTR}/${UQSTR}/historical/${ENSEMBLEMEMBER}/6hrPlevPt/psl/${GRIDDATE}/psl/
  PATHTOFILES2=/glade/collections/cmip/CMIP6/CMIP/${PARENTSTR}/${UQSTR}/historical/${ENSEMBLEMEMBER}/3hr/pr/${GRIDDATE}/pr/
  PATHTOFILES3=/glade/scratch/zarzycki/CMIPTMP/${UQSTR}/
  TOPOFILE=/glade/work/zarzycki/topo-files/CMIP6/${UQSTR}.topo.nc
  THISSED="sed"
elif [[ $(hostname -s) = MET-MAC* ]]; then
  TEMPESTEXTREMESDIR=/Users/cmz5202/Software/tempestextremes/
  PATHTOFILES=/Users/cmz5202/NetCDF/CMIP6/
  TOPOFILE=/Users/cmz5202/NetCDF/topo_files/${UQSTR}.topo.nc
  THISSED="gsed"
elif [[ $(hostname -s) = aci* ]]; then
  TEMPESTEXTREMESDIR=/storage/home/cmz5202/sw/tempestextremes/
  PATHTOFILES=/gpfs/group/cmz5202/default/mjg6459/CMIP6/${PARENTSTR}/${UQSTR}/psl/
  PATHTOFILES3=/storage/home/${LOGNAME}/scratch/CMIPTMP/${UQSTR}/
  TOPOFILE=/gpfs/group/cmz5202/default/topo/${UQSTR}.topo.nc
  TOPOORIGDIR=/storage/home/cmz5202/work/cam_tools/hires-topo/
  THISSED="sed"
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
mkdir -p $(dirname "${TOPOFILE}")

############ TRACK ETCs #####################

## Build trajectory files; the most time intensive part of this.
if [ "$DO_TRACKS" = true ] ; then

  mkdir -p $OUTTRAJDIR

  touch $FILELISTNAME

  #### Files
  find ${PATHTOFILES} -name "psl_6hrPlevPt_${UQSTR}_historical_*_*0101*nc" | sort -n >> $FILELISTNAME
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

  DCU_PSLNAME="psl"    # Name of PSL on netcdf files
  DCU_PSLFOMAG=200.0   # Pressure "anomaly" required for PSL min to be kept
  DCU_PSLFODIST=6.0    # Pressure "anomaly" contour must lie within this distance
  DCU_MERGEDIST=6.0    # Merge two competing PSL minima within this radius
  SN_TRAERA5NGE=6.0     # Maximum distance (GC degrees) ETC can move in successive steps
  SN_TRAJMINLENGTH=10  # Min length of trajectory (nsteps)
  SN_TRAJMAXGAP=3      # Max gap within a trajectory (nsteps)

  STRDETECT="--verbosity 0 --timestride 1 ${CONNECTFLAG} --out cyclones_tempest.${DATESTRING} --closedcontourcmd ${DCU_PSLNAME},${DCU_PSLFOMAG},${DCU_PSLFODIST},0 --mergedist ${DCU_MERGEDIST} --searchbymin ${DCU_PSLNAME} --outputcmd ${DCU_PSLNAME},min,0;PHIS,max,0"
  echo $STRDETECT
  touch cyclones.${DATESTRING}
  ${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STRDETECT} </dev/null
  cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
  rm cyclones_tempest.${DATESTRING}*

  ### Special filtering by time to account for fact that some models have longer PSL datasets than PRECT or PRECSN
  ## Create a temp cyclones file only including dates from STYR onwards, then overwrite full cyclones file before SN
  if [ $STYR \> 0 ]; then
    ${THISSED} -n '/^'${STYR}'/,$p' cyclones.${DATESTRING} > cycfilter.tmp.${DATESTRING}.txt;
    mv cycfilter.tmp.${DATESTRING}.txt cyclones.${DATESTRING}
  fi;

  # Stitch candidate cyclones together
  STRSTITCH="--range ${SN_TRAERA5NGE} --mintime 60h --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist 12.0 --threshold lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"
  ${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_fmt "lon,lat,slp,phis" ${STRSTITCH}

  # Clean up leftover files
  rm cyclones.${DATESTRING}   #Delete candidate cyclone file
  rm ${FILELISTNAME}
  rm log*.txt

  ${TEMPESTEXTREMESDIR}/bin/HistogramNodes --in ${TRAJFILENAME} --nlat 36 --nlon 72 --iloncol 3 --ilatcol 4 --out "${PATHTOFILES3}/density_${ENSEMBLEMEMBER}_${UQSTR}.nc"
fi

############ FIELD EXTRACTION #####################
######## OK, do some node work now that we've tracked

if [ "$DO_EXTRACT" = true ] ; then

  touch $NODELISTNAME

  find ${PATHTOFILES2} -name "pr_3hr_${UQSTR}_historical_${ENSEMBLEMEMBER}*.nc" | sort -n >> $NODELISTNAME

  ## Create files and folder to hold filtered data
  cp $NODELISTNAME $NODELISTNAMEFILT
  $THISSED -i 's/pr_/pr_FILT_/g' $NODELISTNAMEFILT
  $THISSED -i 's?'${PATHTOFILES2//\?/\.}'?'${PATHTOFILES3}'?g' $NODELISTNAMEFILT
  mkdir -p ${PATHTOFILES3}

  ${TEMPESTEXTREMESDIR}/bin/NodeFileFilter --in_nodefile ${TRAJFILENAME} --in_fmt "lon,lat,slp,phis" --in_data_list ${NODELISTNAME} --out_data_list ${NODELISTNAMEFILT} --var "pr" --bydist 25.0 --maskvar "mask" #  --nearbyblobs "pr,2.0,>=,5.0e-5,50.0"     pr is kg/m2/s so pr/1000 => m/s

#    TASFILES=`find ${PATHTOFILES2//pr/tas} -name "tas_3hr_${UQSTR}_historical_${ENSEMBLEMEMBER}*.nc" | sort -n`
#    for f in $TASFILES; do
#      ## Find matching *FILT* file with last 15 characters (date and nc)
#      THISPRECFILE=`ls ${PATHTOFILES3}/*FILT*${f:(-15)}`
#      echo "Appending... $f to $THISPRECFILE"
#      ## Append snow to FILT file
#      ncks -A -v tas $f $THISPRECFILE
#      ## Use mask var to mask off non-cyclone
#      ncap2 -O -s 'where(mask!=1) tas=0.0' $THISPRECFILE $THISPRECFILE
#    done

   SNOWFILES=`find ${PATHTOFILES2//pr/prsn} -name "prsn_3hr_${UQSTR}_historical_${ENSEMBLEMEMBER}*.nc" | sort -n`
   for f in $SNOWFILES; do
     ## Find matching *FILT* file with last 15 characters (date and nc)
     THISPRECFILE=`ls ${PATHTOFILES3}/*FILT*${f:(-15)}`
     echo "Appending... $f to $THISPRECFILE"
     ## Append snow to FILT file
     ncks -A -v prsn $f $THISPRECFILE
     ## Use mask var to mask off non-cyclone
     ncap2 -O -s 'where(mask!=1) prsn=0.0' $THISPRECFILE $THISPRECFILE
   done

  pattern=${PATHTOFILES3}"/pr_FILT_3hr_${UQSTR}_historical_${ENSEMBLEMEMBER}*nc"
  dumfiles=( $pattern )
  ${TEMPESTEXTREMESDIR}/bin/NodeFileEditor --in_nodefile ${TRAJFILENAME} --in_data ${dumfiles[0]} --in_fmt "lon,lat,slp,phis" --out_nodefile ${TRAJFILENAME}_strong.txt --out_fmt "lon,lat,slp,phis" --colfilter "slp,<=,99000."

### KEEP COMMENTED
# in hypothetical world, where flatten time would go
#  for f in /glade/u/home/zarzycki/scratch/CMIPTMP/NESM3/pr_*nc
#  do
#    echo $f
#    ncl flatten-time.ncl 'filename="'${f}'"'
#  done

  ${TEMPESTEXTREMESDIR}/bin/NodeFileCompose --in_nodefile ${TRAJFILENAME}_strong.txt --in_fmt "lon,lat,slp,phis" --in_data_list "${NODELISTNAMEFILT}" --var "pr" --max_time_delta "2h" --out_data "${PATHTOFILES3}/composite_${ENSEMBLEMEMBER}_${UQSTR}.nc" --dx 1.0 --resx 80 --op "mean" --snapshots
  
  # Clean up
  rm ${NODELISTNAME}
  rm ${NODELISTNAMEFILT}

fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime},${TRAJFILENAME}\n" >> timing.txt
