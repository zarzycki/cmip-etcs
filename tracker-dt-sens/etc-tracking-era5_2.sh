#!/bin/bash -l

############ USER OPTIONS #####################

## If using unstructured data, need connect file, otherwise empty
CONNECTFLAG=""
DO_TRACKS=true
DO_EXTRACT=false
ENSEMBLEMEMBER=r1i1p1f1   #r1, r2, etc.
UQSTR="ERA5"

## Unique string

#UQSTR="GFDL-CM4"
#PARENTSTR="NOAA-GFDL"
#GRIDDATE="gr1/v20180701"
#STYR="-1"

#UQSTR="GISS-E2-1-G"
#PARENTSTR="NASA-GISS"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="CMCC-CM2-HR4"
#PARENTSTR="CMCC"
#GRIDDATE="gn/v????????"
#STYR="-1"

### CMZ: issues with snow files
#UQSTR="IPSL-CM6A-LR"
#PARENTSTR="IPSL"
#GRIDDATE="gr/v????????"
#STYR="-1"

#UQSTR="MPI-ESM-1-2-HAM"
#PARENTSTR="HAMMOZ-Consortium"
#GRIDDATE="gr/v????????"
#STYR="-1"

#UQSTR="KIOST-ESM"
#PARENTSTR="KIOST"
#GRIDDATE="gr/v????????"
#STYR="-1"

#UQSTR="KACE-1-0-G"
#PARENTSTR="NIMS-KMA"
#GRIDDATE="gr/v????????"
#STYR="-1"

#UQSTR="EC-Earth3"
#PARENTSTR="EC-Earth-Consortium"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="AWI-ESM-1-1-LR"
#PARENTSTR="AWI"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="MIROC6"
#PARENTSTR="MIROC"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="MRI-ESM2-0"
#PARENTSTR="MRI"
#GRIDDATE="gn/v????????"
#STYR="-1"

## CMZ: issues with time coordinate in pr/prsn files
#UQSTR="BCC-CSM2-MR"
#PARENTSTR="BCC"
#GRIDDATE="gn/v????????"
#STYR="-1"

# CMZ: need to flatten time
#UQSTR="NESM3"
#PARENTSTR="NUIST"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="SAM0-UNICON"
#PARENTSTR="SNU"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="NorESM2-LM"
#PARENTSTR="NCC"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="TaiESM1"
#PARENTSTR="AS-RCEC"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="CMCC-CM2-SR5"
#PARENTSTR="CMCC"
#GRIDDATE="gn/v????????"
#STYR="-1"

#UQSTR="ACCESS-ESM1-5"
#PARENTSTR="CSIRO"
#GRIDDATE="gn/v????????"
#STYR="-1"

############ MACHINE SPECIFIC AUTO-CONFIG #####################

#if [[ $HOSTNAME = casper* ]]; then
  TEMPESTEXTREMESDIR=/glade/u/home/zarzycki/work/tempestextremes_noMPI/
  PATHTOFILES=/glade/derecho/scratch/mgore/ERA_data/
  PATHTOFILES2=/glade/campaign/collections/rda/data/ds633.0/e5.oper.fc.sfc.meanflux/198001/
  PATHTOFILES3=/glade/derecho/scratch/mgore/CMIPTMP/
  TOPOORIGDIR=/glade/u/home/zarzycki/work/cam_tools/hires-topo/
  TOPOFILE=/glade/work/mgore/topo-files/ERA5.topo.nc
  THISSED="sed"
#elif [[ $(hostname -s) = MET-MAC* ]]; then
#  TEMPESTEXTREMESDIR=/Users/cmz5202/Software/tempestextremes/
#  PATHTOFILES=/Users/cmz5202/NetCDF/CMIP6/
#  TOPOFILE=/Users/cmz5202/NetCDF/topo_files/${UQSTR}.topo.nc
#  THISSED="gsed"
#elif [[ $(hostname -s) = aci* ]]; then
#  TEMPESTEXTREMESDIR=/storage/home/cmz5202/sw/tempestextremes/
#  PATHTOFILES=/gpfs/group/cmz5202/default/mjg6459/HighResMIP/${PARENTSTR}/${UQSTR}/
#  PATHTOFILES3=/storage/home/${LOGNAME}/scratch/CMIPTMP/${UQSTR}/
#  TOPOFILE=/storage/home/mjg6459/work/topo/${UQSTR}.topo.nc
##  TOPOFILE=/gpfs/group/cmz5202/default/topo/${UQSTR}.topo.nc
#  TOPOORIGDIR=/storage/home/cmz5202/work/cam_tools/hires-topo/
#  THISSED="sed"
#elif [[ $(hostname -s) = submit* ]]; then
#  TEMPESTEXTREMESDIR=/storage/home/cmz5202/sw/tempestextremes/
#  PATHTOFILES=/gpfs/group/cmz5202/default/mjg6459/HighResMIP/${PARENTSTR}/${UQSTR}/
#  PATHTOFILES3=/storage/home/${LOGNAME}/scratch/CMIPTMP/${UQSTR}/
#  TOPOFILE=/storage/home/mjg6459/work/topo/${UQSTR}.topo.nc
##  TOPOFILE=/gpfs/group/cmz5202/default/topo/${UQSTR}.topo.nc
#  TOPOORIGDIR=/storage/home/cmz5202/work/cam_tools/hires-topo/
#  THISSED="sed"
#else
#  echo "Can't figure out hostname, exiting"
#  exit
#fi

############ TRACKER MECHANICS #####################
starttime=$(date -u +"%s")

DATESTRING=`date +"%s%N"`
OUTTRAJDIR="/glade/derecho/scratch/zarzycki/ETC-track/trajs/"
TRAJFILENAME=${OUTTRAJDIR}/trajsens2.txt
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
  find ${PATHTOFILES} -name "NorthAm_ERA5_MSLP_1980.nc" | sort -n >> $FILELISTNAME
#  find ${PATHTOFILES} -name "E3SMHRHR.h1.*nc" | sort -n >> $FILELISTNAME
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

#  Original
#  DCU_PSLNAME="psl"    # Name of PSL on netcdf files
#  DCU_PSLFOMAG=200.0   # Pressure "anomaly" required for PSL min to be kept
#  DCU_PSLFODIST=6.0    # Pressure "anomaly" contour must lie within this distance
#  DCU_MERGEDIST=6.0    # Merge two competing PSL minima within this radius
#  SN_TRAERA5NGE=6.0     # Maximum distance (GC degrees) ETC can move in successive steps
#  SN_TRAJMINLENGTH=10  # Min length of trajectory (nsteps)
#  SN_TRAJMAXGAP=3      # Max gap within a trajectory (nsteps)

#  Sensitivity 1
#  DCU_PSLNAME="MSL"    # Name of PSL on netcdf files
#  DCU_PSLFOMAG=200.0   # Pressure "anomaly" required for PSL min to be kept
#  DCU_PSLFODIST=6.0    # Pressure "anomaly" contour must lie within this distance
#  DCU_MERGEDIST=6.0    # Merge two competing PSL minima within this radius
#  SN_TRAJRANGE=6.0     # Maximum distance (GC degrees) ETC can move in successive steps
#  SN_TRAJMINLENGTH=4   # Min length of trajectory (nsteps)
#  SN_TRAJMAXGAP=1      # Max gap within a trajectory (nsteps)
#  SN_MINENDPOINT=9.01  # Min travel distance
#  SN_MINTIME="24h"

#  Sensitivity 2
  DCU_PSLNAME="MSL"    # Name of PSL on netcdf files
  DCU_PSLFOMAG=200.0   # Pressure "anomaly" required for PSL min to be kept
  DCU_PSLFODIST=6.0    # Pressure "anomaly" contour must lie within this distance
  DCU_MERGEDIST=9.01   # Merge two competing PSL minima within this radius
  SN_TRAJRANGE=1.5   # Maximum distance (GC degrees) ETC can move in successive steps
  SN_TRAJMINLENGTH=4  # Min length of trajectory (nsteps)
  SN_TRAJMAXGAP=1      # Max gap within a trajectory (nsteps)
  SN_MINENDPOINT=12.0  # Min travel distance
  SN_MINTIME="24h"

  STRDETECT="--verbosity 0 --timestride 6 ${CONNECTFLAG} --out cyclones_tempest.${DATESTRING} --closedcontourcmd ${DCU_PSLNAME},${DCU_PSLFOMAG},${DCU_PSLFODIST},0 --mergedist ${DCU_MERGEDIST} --searchbymin ${DCU_PSLNAME} --outputcmd ${DCU_PSLNAME},min,0;PHIS,max,0"
  echo $STRDETECT
  touch cyclones.${DATESTRING}
  ${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STRDETECT} </dev/null
  cat cyclones_tempest.* >> cyclones.${DATESTRING}
  rm cyclones_tempest.${DATESTRING}*

  ### Special filtering by time to account for fact that some models have longer PSL datasets than PRECT or PRECSN
  ## Create a temp cyclones file only including dates from STYR onwards, then overwrite full cyclones file before SN
  if [ $STYR \> 0 ]; then
    ${THISSED} -n '/^'${STYR}'/,$p' cyclones.${DATESTRING} > cycfilter.tmp.${DATESTRING}.txt;
    mv cycfilter.tmp.${DATESTRING}.txt cyclones.${DATESTRING}
  fi;

  exit

  # Stitch candidate cyclones together
#  STRSTITCH="--range ${SN_TRAERA5NGE} --mintime 60h --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist 12.0 --threshold lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"
  STRSTITCH="--range ${SN_TRAJRANGE} --mintime ${SN_MINTIME} --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist ${SN_MINENDPOINT} --threshold lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"
  ${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_fmt "lon,lat,psl,phis" ${STRSTITCH}

  # Clean up leftover files
  #rm cyclones.${DATESTRING}   #Delete candidate cyclone file
  rm ${FILELISTNAME}
  #rm log*.txt

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

  ${TEMPESTEXTREMESDIR}/bin/NodeFileFilter --in_nodefile ${TRAJFILENAME} --in_fmt "lon,lat,psl,phis" --in_data_list ${NODELISTNAME} --out_data_list ${NODELISTNAMEFILT} --var "pr" --bydist 25.0 --maskvar "mask" #  --nearbyblobs "pr,2.0,>=,5.0e-5,50.0"     pr is kg/m2/s so pr/1000 => m/s

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
  ${TEMPESTEXTREMESDIR}/bin/NodeFileEditor --in_nodefile ${TRAJFILENAME} --in_data ${dumfiles[0]} --in_fmt "lon,lat,psl,phis" --out_nodefile ${TRAJFILENAME}_strong.txt --out_fmt "lon,lat,slp,phis" --colfilter "slp,<=,99000."

### KEEP COMMENTED
# in hypothetical world, where flatten time would go
#  for f in /glade/u/home/zarzycki/scratch/CMIPTMP/NESM3/pr_*nc
#  do
#    echo $f
#    ncl flatten-time.ncl 'filename="'${f}'"'
#  done

  ${TEMPESTEXTREMESDIR}/bin/NodeFileCompose --in_nodefile ${TRAJFILENAME}_strong.txt --in_fmt "lon,lat,psl,phis" --in_data_list "${NODELISTNAMEFILT}" --var "pr" --max_time_delta "2h" --out_data "${PATHTOFILES3}/composite_${ENSEMBLEMEMBER}_${UQSTR}.nc" --dx 1.0 --resx 80 --op "mean" --snapshots

  # Clean up
  rm ${NODELISTNAME}
  rm ${NODELISTNAMEFILT}

fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime},${TRAJFILENAME}\n" >> timing.txt
