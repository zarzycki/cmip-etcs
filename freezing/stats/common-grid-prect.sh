#!/bin/bash

THISDIR=$PWD
FZRADIR="/glade/derecho/scratch/zarzycki/FZRA/"

# Define years to process
START_YEAR=1980
END_YEAR=2020

cd $FZRADIR/precip/ || { echo "Failed to change directory to $FZRADIR/precip/. Exiting."; exit 1; }
pwd

# Define datasets and their input file patterns
declare -A datasets=(
    ["CFSR"]="CFSR/CFSR.PREC.YYYY.nc"
    ["ERA5"]="ERA5/NorthAm_ERAML_PREC_YYYY.nc"
    ["CR20"]="CR20/anl_mean_YYYY_PR.nc"
    ["JRA"]="JRA/fcst_phy2m125_tprat_YYYY.nc"
    ["IMERG"]="IMERG/3IMERG.6hr.YYYY.nc"
)

# Loop through each dataset
for dataset in "${!datasets[@]}"; do
    echo "Processing ${dataset} dataset..."
    pattern="${datasets[$dataset]}"

    # Loop through years
    for year in $(seq $START_YEAR $END_YEAR); do
        input_file="${pattern/YYYY/$year}"
        output_file="pre-${dataset}-${year}.nc"
        echo "Processing ${dataset} data for year ${year}..."

        # Check if input file exists
        if [ ! -f "$input_file" ]; then
            echo "Warning: Input file not found: ${input_file}"
            continue
        fi

        # Dataset-specific processing
        case $dataset in
            "CFSR")
                ncremap -a patch -d $THISDIR/target.nc -i "$input_file" -o "$output_file"
                ncrename -O -v PRATE,PRECT "$output_file" "$output_file"
                ncap2 -O -s 'PRECT=PRECT*86400' -s 'PRECT@units="mm/day"' "$output_file" "$output_file"
                ;;

            "ERA5")
                ncremap -a patch -d $THISDIR/target.nc -i "$input_file" -o "$output_file"
                ncks -O -d time,,,6 "$output_file" "$output_file"
                ncrename -O -v MTPR,PRECT "$output_file" "$output_file"
                ncap2 -O -s 'PRECT=PRECT*86400' -s 'PRECT@units="mm/day"' "$output_file" "$output_file"
                ;;

            "CR20")
                ncremap -a patch -d $THISDIR/target.nc -i "$input_file" -o "$output_file"
                ncks -O -3 "$output_file" "$output_file"
                ncrename -O -d .latitude,lat -d .longitude,lon "$output_file" "$output_file"
                ncrename -O -v prate,PRECT "$output_file" "$output_file"
                ncap2 -O -s 'PRECT=PRECT*86400' -s 'PRECT@units="mm/day"' "$output_file" "$output_file"
                #ncks -O -4 -L 1 "$output_file" "$output_file"
                ;;

            "JRA")
                # Create temporary file for JRA processing
                ncks -d forecast_time1,0,0 "$input_file" tmp.nc
                ncwa -O -a forecast_time1 tmp.nc tmp.nc
                ncks -O -3 tmp.nc tmp.nc
                ncrename -O -d .initial_time0_hours,time -v .initial_time0_hours,time \
                        -d .g0_lat_2,lat -d .g0_lon_3,lon \
                        -v .g0_lat_2,lat -v .g0_lon_3,lon tmp.nc tmp.nc
                ncremap -a patch -d $THISDIR/target.nc -i tmp.nc -o "$output_file"
                ncrename -O -v TPRAT_GDS0_SFC_ave3h,PRECT "$output_file" "$output_file"
                #ncks -O -4 -L 1 "$output_file" "$output_file"
                ncatted -O -a units,PRECT,m,c,"mm/day" "$output_file" "$output_file"
                rm -fv tmp.nc
                ;;

            "IMERG")
                ncremap -a patch -d $THISDIR/target.nc -i "$input_file" -o "$output_file"
                # Convert from mm/6hr to mm/day
                ncap2 -O -s "prec=4*prec" "$output_file" "$output_file"
                ncrename -O -v prec,PRECT "$output_file" "$output_file"
                ncatted -O -a units,PRECT,m,c,"mm/day" "$output_file" "$output_file"
                ;;
        esac

        if [ $? -eq 0 ]; then
            echo "Successfully processed ${dataset} ${year}"
        else
            echo "Error processing ${dataset} ${year}"
        fi
    done
done
