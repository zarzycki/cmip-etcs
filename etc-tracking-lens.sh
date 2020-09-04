#!/bin/bash -l

############ USER OPTIONS #####################

## If using unstructured data, need connect file, otherwise empty
CONNECTFLAG="" 
DO_TRACKS=false
DO_EXTRACT=true

## Unique string
UQSTR="CESM1-LENS"
PARENTSTR="CESM"

############ MACHINE SPECIFIC AUTO-CONFIG #####################

if [[ $HOSTNAME = cheyenne* ]]; then
  TEMPESTEXTREMESDIR=/glade/u/home/zarzycki/work/tempestextremes_noMPI/
  PATHTOFILES=/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/PSL/
  PATHTOFILES2=/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/PRECT/
  PATHTOFILES3=/glade/scratch/zarzycki/CMIPTMP/${UQSTR}/
  TOPOFILE=/glade/work/zarzycki/topo-files/CMIP6/${UQSTR}.topo.nc
  THISSED="sed"
elif [[ $(hostname -s) = MET-MAC* ]]; then
  TEMPESTEXTREMESDIR=/Users/cmz5202/Software/tempestextremes/
  PATHTOFILES=/Users/cmz5202/NetCDF/CMIP6/
  TOPOFILE=/Users/cmz5202/NetCDF/topo_files/${UQSTR}.topo.nc
  THISSED="gsed"
else
  echo "Can't figure out hostname, exiting"
  exit
fi

############ TRACKER MECHANICS #####################
starttime=$(date -u +"%s")

DATESTRING=`date +"%s%N"`
OUTTRAJDIR="./trajs/"
TRAJFILENAME=${OUTTRAJDIR}/trajectories.txt.${UQSTR}
FILELISTNAME=filelist.txt.${DATESTRING}
NODELISTNAME=nodelist.txt.${DATESTRING}
NODELISTNAMEFILT=nodelistfilt.txt.${DATESTRING}

####

mkdir -p ${PATHTOFILES3}

############ TRACK ETCs #####################

## Build trajectory files; the most time intensive part of this.
if [ "$DO_TRACKS" = true ] ; then

  mkdir -p $OUTTRAJDIR

  touch $FILELISTNAME

  #### Files
  find ${PATHTOFILES} -name "b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PSL.1990010100Z-2005123118Z.nc" | sort -n >> $FILELISTNAME
  # Add any static file(s) to each line
  $THISSED -e 's?$?;'"${TOPOFILE}"'?' -i $FILELISTNAME
  cat $FILELISTNAME

  # Check if topo file exists, if not generate one using the first PSL file
  if [ ! -f ${TOPOFILE} ]; then
    echo "Topo file ${TOPOFILE} not found!"
    FIRSTPSLFILE=$(echo $(head -n 1 ${FILELISTNAME}) | awk -F';' '{print $1}')
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
  # used for TCs, not ETCs (right now)
  #SN_MAXTOPO=150.0   
  #SN_MINWIND=10.0

  STRDETECT="--verbosity 0 --timestride 1 ${CONNECTFLAG} --out cyclones_tempest.${DATESTRING} --closedcontourcmd ${DCU_PSLNAME},${DCU_PSLFOMAG},${DCU_PSLFODIST},0 --mergedist ${DCU_MERGEDIST} --searchbymin ${DCU_PSLNAME} --outputcmd ${DCU_PSLNAME},min,0;${DCU_PSLNAME},min,0;PHIS,max,0"
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
  STRSTITCH="--range ${SN_TRAJRANGE} --mintime 60h --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist 12.0 --threshold lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"
  ${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_fmt "lon,lat,slp,slp,phis" ${STRSTITCH}

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

  find ${PATHTOFILES2} -name "b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PRECT.1990010100Z-2005123118Z.nc" | sort -n >> $NODELISTNAME

  ## Create files and folder to hold filtered data
  cp $NODELISTNAME $NODELISTNAMEFILT
  $THISSED -i 's?'${PATHTOFILES2//\?/\.}'?'${PATHTOFILES3}'?g' $NODELISTNAMEFILT
  $THISSED -i 's/PRECT./PRECT_FILT./g' $NODELISTNAMEFILT
  mkdir -p ${PATHTOFILES3}
  
  
  #${TEMPESTEXTREMESDIR}/bin/NodeFileFilter --in_nodefile ${TRAJFILENAME} --in_fmt "lon,lat,slp,slp,phis" --in_data_list ${NODELISTNAME} --out_data_list ${NODELISTNAMEFILT} --var "PRECT" --bydist 25.0 --maskvar "mask" #  --nearbyblobs "pr,2.0,>=,5.0e-5,50.0"     pr is kg/m2/s so pr/1000 => m/s

  # Declare an array of string with type
  declare -a StringArray=("TREFHT" "CLDTOT" "QREFHT" "TGCLDIWP" "TGCLDLWP")

  EXTRASTR=""
  # Iterate the string array using for loop
  for val in ${StringArray[@]}; do
    VAR=$val
    EXTRASTR=$EXTRASTR","$VAR
    SNOWFILES=`find ${PATHTOFILES2//PRECT/${VAR}} -name "b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.${VAR}.1990010100Z-2005123118Z.nc" | sort -n`
    echo $SNOWFILES
    for f in $SNOWFILES; do
      echo "Appending... $f"
      # Find matching *FILT* file with last 15 characters (date and nc)
      THISPRECFILE=`ls ${PATHTOFILES3}/*FILT*${f:(-15)}`
      # Append snow to FILT file
      ##ncks -A -v ${VAR} $f $THISPRECFILE
      # Use mask var to mask off non-cyclone
      echo "... masking"
      ##ncap2 -O -s 'where(mask!=1) '${VAR}'=0.0' $THISPRECFILE $THISPRECFILE
    done
  done

  pattern=${PATHTOFILES3}"/b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PRECT_FILT.1990010100Z-2005123118Z.nc"
  dumfiles=( $pattern )
  ${TEMPESTEXTREMESDIR}/bin/NodeFileEditor --in_nodefile ${TRAJFILENAME} --in_data ${dumfiles[0]} --in_fmt "lon,lat,slp,slp,phis" --out_nodefile ${TRAJFILENAME}_strong.txt --out_fmt "lon,lat,slp,slp,phis" --colfilter "slp,<=,99000."

  ${TEMPESTEXTREMESDIR}/bin/NodeFileCompose --in_nodefile ${TRAJFILENAME}_strong.txt --in_fmt "lon,lat,slp,slp,phis" --in_data_list "${NODELISTNAMEFILT}" --var "PRECT"${EXTRASTR} --max_time_delta "2h" --out_data "${PATHTOFILES3}/composite_${UQSTR}.nc" --dx 1.0 --resx 80 --op "mean" --snapshots
  
  # Clean up
  rm ${NODELISTNAME}
  rm ${NODELISTNAMEFILT}

fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime},${TRAJFILENAME}\n" >> timing.txt


