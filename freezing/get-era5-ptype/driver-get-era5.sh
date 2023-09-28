#!/bin/bash

# NOTE, need cdsapi Python
#(base) zarzycki@cheyenne3:~/scratch/ptype> module load conda
#(base) zarzycki@cheyenne3:~/scratch/ptype> conda activate npl
# conda activate meteo

## declare an array variable
declare -a months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

STYR=1980
ENYR=2022
OUTDIR=/glade/u/home/zarzycki/scratch/ptype/

mkdir -p -v $OUTDIR

## now loop through the above array
for jj in $(seq $STYR $ENYR);
do
  for ii in "${months[@]}"
  do
    echo $jj' '$ii
    OUTFILE=${OUTDIR}/ptype.${jj}.${ii}.nc
    if [ -f "$OUTFILE" ]; then
      echo "$OUTFILE exists, not downloading"
    else
      echo "$OUTFILE doesn't exist..."
      python get-data.py "${jj}" "${ii}" "${OUTDIR}"
      echo "Compressing..."
      ncks -O -4 -L 1 ${OUTFILE} ${OUTFILE}
    fi
  done
done
