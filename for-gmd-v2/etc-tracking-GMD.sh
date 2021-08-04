#!/bin/bash -l

############ USER OPTIONS #####################

## If using unstructured data, need connect file, otherwise empty
CONNECTFLAG="" 
DO_TRACKS=false
DO_EXTRACT=true

## Unique string
UQSTR="CESM2-LENS"
PARENTSTR="CESM"
FILTVAR=QREFHT
############ MACHINE SPECIFIC AUTO-CONFIG #####################

if [[ $HOSTNAME = cheyenne* ]]; then
  TEMPESTEXTREMESDIR=/glade/u/home/zarzycki/work/tempestextremes_noMPI/
  PATHTOFILES=/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/PSL/
  PATHTOFILES2=/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/${FILTVAR}/
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
FILELISTNAME=B20TRC5CNBDRD.001.PS_list.txt
NODELISTNAME=B20TRC5CNBDRD.001.${FILTVAR}_list.txt
NODELISTNAMEFILT=B20TRC5CNBDRD.001.${FILTVAR}_FILT_list.txt

####

mkdir -p ${PATHTOFILES3}

############ TRACK ETCs #####################

## Build trajectory files; the most time intensive part of this.
if [ "$DO_TRACKS" = true ] ; then

  mkdir -p $OUTTRAJDIR
  
  rm -v B20TRC5CNBDRD.001.PS_list.txt

  #### Files
  find ${PATHTOFILES} -name "b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.PSL.1990010100Z-2005123118Z.nc" | sort -n >> B20TRC5CNBDRD.001.PS_list.txt

  ${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data_list B20TRC5CNBDRD.001.PS_list.txt --out cyclone_candidates --closedcontourcmd PSL,200.0,6.0,0 --mergedist 6.0 --searchbymin "PSL" --outputcmd "PSL,min,0"
  
  ls cyclone_candidates* > candidate_list.txt

  ${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_fmt "lon,lat,slp" --caltype "noleap" --in_list candidate_list.txt --out etc-all-traj.txt --range 6.0 --mintime 60h --maxgap 18h --min_endpoint_dist 12.0 --threshold "lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"
  
fi

############ FIELD EXTRACTION #####################
######## OK, do some node work now that we've tracked

if [ "$DO_EXTRACT" = true ] ; then

  rm $NODELISTNAME $NODELISTNAMEFILT

  find ${PATHTOFILES2} -name "b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h2.${FILTVAR}.1990010100Z-2005123118Z.nc" | sort -n >> $NODELISTNAME

  ## Create files and folder to hold filtered data
  cp $NODELISTNAME $NODELISTNAMEFILT
  $THISSED -i 's?'${PATHTOFILES2//\?/\.}'?'${PATHTOFILES3}'?g' $NODELISTNAMEFILT
  $THISSED -i 's/${FILTVAR}./${FILTVAR}_FILT./g' $NODELISTNAMEFILT
  mkdir -p ${PATHTOFILES3}
  
${TEMPESTEXTREMESDIR}/bin/NodeFileFilter --in_nodefile etc-all-traj.txt --in_fmt "lon,lat,slp" --in_data_list B20TRC5CNBDRD.001.${FILTVAR}_list.txt --out_data_list B20TRC5CNBDRD.001.${FILTVAR}_FILT_list.txt --var "${FILTVAR}" --bydist 25.0 --maskvar "mask"
  
${TEMPESTEXTREMESDIR}/bin/NodeFileEditor --in_nodefile etc-all-traj.txt --in_data_list B20TRC5CNBDRD.001.${FILTVAR}_list.txt --in_fmt "lon,lat,slp" --out_fmt "lon,lat,slp" --out_nodefile etc-strong-traj.txt --colfilter "slp,<=,99000."
  
${TEMPESTEXTREMESDIR}/bin/NodeFileCompose --in_nodefile etc-strong-traj.txt --in_fmt "lon,lat,slp" --in_data_list B20TRC5CNBDRD.001.${FILTVAR}_FILT_list.txt --out_data "composite_${FILTVAR}.nc" --var "${FILTVAR}" --max_time_delta "2h" --dx 1.0 --resx 80 --op "mean" #--snapshots 

fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime},${TRAJFILENAME}\n" >> timing.txt


