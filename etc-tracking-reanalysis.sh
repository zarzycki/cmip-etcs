#!/bin/bash -l

############ USER OPTIONS #####################

## Unique string
UQSTR=GISS-E2-1-G

## If using unstructured data, need connect file, otherwise empty
CONNECTFLAG="" 

############ MACHINE SPECIFIC AUTO-CONFIG #####################

if [[ $(hostname -s) = cheyenne* ]]; then
  TEMPESTEXTREMESDIR=/Users/cmz5202/Software/tempestextremes/
  PATHTOFILES=/Users/cmz5202/NetCDF/CMIP6/
  TOPOFILE=/Users/cmz5202/NetCDF/topo_files/${UQSTR}.topo.nc
elif [[ $(hostname -s) = MET-MAC* ]]; then
  TEMPESTEXTREMESDIR=/Users/cmz5202/Software/tempestextremes/
  PATHTOFILES=/Users/cmz5202/NetCDF/CMIP6/
  TOPOFILE=/Users/cmz5202/NetCDF/topo_files/${UQSTR}.topo.nc
else
  echo "Can't figure out hostname, exiting"
  exit
fi

############ TRACKER MECHANICS #####################
starttime=$(date -u +"%s")

DATESTRING=`date +"%s%N"`
FILELISTNAME=filelist.txt.${DATESTRING}
OUTTRAJDIR="./trajs/"
TRAJFILENAME=${OUTTRAJDIR}/trajectories.txt.${UQSTR}

mkdir -p $OUTTRAJDIR

touch $FILELISTNAME

#for zz in ${YEARSTRARR}
#do
  find ${PATHTOFILES} -name "psl_6hrPlevPt_${UQSTR}_historical_*_*0101*nc" | grep -v 1950 | sort -n >> $FILELISTNAME
#done
# Add any static file(s) to each line
gsed -e 's?$?;'"${TOPOFILE}"'?' -i $FILELISTNAME

DCU_PSLNAME="psl"    # Name of PSL on netcdf files
DCU_PSLFOMAG=200.0   # Pressure "anomaly" required for PSL min to be kept
DCU_PSLFODIST=6.0    # Pressure "anomaly" contour must lie within this distance
DCU_MERGEDIST=6.0    # Merge two competing PSL minima within this radius
SN_TRAERA5NGE=6.0     # Maximum distance (GC degrees) ETC can move in successive steps
SN_TRAJMINLENGTH=10  # Min length of trajectory (nsteps)
SN_TRAJMAXGAP=3      # Max gap within a trajectory (nsteps)
# used for TCs, not ETCs (right now)
#SN_MAXTOPO=150.0   
#SN_MINWIND=10.0

# Detect candidates
STRDETECT="--verbosity 0 --timestride 1 ${CONNECTFLAG} --out cyclones_tempest.${DATESTRING} --closedcontourcmd ${DCU_PSLNAME},${DCU_PSLFOMAG},${DCU_PSLFODIST},0 --mergedist ${DCU_MERGEDIST} --searchbymin ${DCU_PSLNAME} --outputcmd ${DCU_PSLNAME},min,0;${DCU_PSLNAME},min,0;PHIS,max,0"
echo $STRDETECT
touch cyclones.${DATESTRING}
${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STRDETECT} </dev/null
cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}

# Stitch candidate cyclones together
STRSTITCH="--range ${SN_TRAERA5NGE} --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --min_endpoint_dist 12.0 --threshold lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"
${TEMPESTEXTREMESDIR}/bin/StitchNodes --format "i,j,lon,lat,slp,slp,phis" ${STRSTITCH}

# Clean up
rm cyclones.${DATESTRING}   #Delete candidate cyclone file
rm cyclones_tempest.${DATESTRING}*
rm ${FILELISTNAME}
rm log*.txt

${TEMPESTEXTREMESDIR}/bin/HistogramNodes --in ${TRAJFILENAME} --nlat 45 --nlon 90 --iloncol 3 --ilatcol 4 --out test.nc

######## OK, do some node work now that we've tracked


NODELISTNAME=nodelist.txt.${DATESTRING}
NODELISTNAMEFILT=nodelistfilt.txt.${DATESTRING}
touch $NODELISTNAME
#for zz in ${YEARSTRARR}
#do
  find ${PATHTOFILES} -name "pr_3hr_${UQSTR}_historical_r1i1p1f1*.nc_CORR.nc" | sort -n >> $NODELISTNAME
#done

cp $NODELISTNAME $NODELISTNAMEFILT
sed -i '' 's/pr_/pr_FILT_/g' $NODELISTNAMEFILT

#######

${TEMPESTEXTREMESDIR}/bin/NodeFileFilter --in_nodefile ${TRAJFILENAME} --in_fmt "lon,lat,slp,slp,phis" --in_data_list ${NODELISTNAME} --out_data_list ${NODELISTNAMEFILT} --var "pr" --nearbyblobs "pr,2.0,>=,10.0e-5,40.0" --maskvar "mask"

pattern="/Users/cmz5202/NetCDF/CMIP6/pr_FILT_3hr_${UQSTR}_historical_r1i1p1f1*nc"
dumfiles=( $pattern )
${TEMPESTEXTREMESDIR}/bin/NodeFileEditor --in_nodefile ${TRAJFILENAME} --in_data ${dumfiles[0]} --in_fmt "lon,lat,slp,slp,phis" --out_nodefile ${TRAJFILENAME}_strong.txt --out_fmt "lon,lat,slp,slp,phis" --col_filter "slp,<=,97000."

${TEMPESTEXTREMESDIR}/bin/NodeFileCompose --in_nodefile ${TRAJFILENAME}_strong.txt --in_fmt "lon,lat,slp,slp,phis" --in_data_list "${NODELISTNAMEFILT}" --var "pr" --out_data composite.nc --dx 1.0 --resx 80 --op "mean" --snapshots

# Clean up
rm ${NODELISTNAME}
rm ${NODELISTNAMEFILT}

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime},${TRAJFILENAME}\n" >> timing.txt


