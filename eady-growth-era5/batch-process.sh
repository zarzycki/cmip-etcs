#!/bin/bash -l

#PBS -N eady_ERA5
#PBS -A P93300642
#PBS -l select=1:ncpus=36:mem=240GB
#PBS -l walltime=12:00:00
#PBS -q casper
#PBS -j oe

module load parallel
module load ncl
module load peak_memusage

# use numcores = 32 for 1deg runs, 16 for high-res runs
NUMCORES=36
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

for YEAR in {1980..2019} ; do

  if [ `expr $YEAR % 400` -eq 0 ] ; then
    SEQ=365
  elif [ `expr $YEAR % 100` -eq 0 ] ; then
    SEQ=364
  elif [ `expr $YEAR % 4` -eq 0 ] ; then
    SEQ=365
  else
    SEQ=364
  fi

  echo $SEQ

  for j in `seq 0 ${SEQ}`   # 0 36
  do
    st=$j
    en=$j
    for i in `seq ${st} ${en}`
    do
      yyyy=`date -d "${YEAR}-01-01 $i days" +%Y`
      mm=`date -d "${YEAR}-01-01 $i days" +%m`
      dd=`date -d "${YEAR}-01-01 $i days" +%d`
      THISDATE=${yyyy}${mm}${dd}
      
      LINECOMMAND="cd /glade/u/home/zarzycki/sw/cmip-etcs/eady-growth-era5/ ;
          ncl -n eady.ncl
          'YYYYMMDD=\"'${THISDATE}'\"'
          "
          
      echo ${LINECOMMAND} >> ${COMMANDFILE}
    done
  done

done

peak_memusage.exe parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

# Cleanup
rm ${COMMANDFILE}
