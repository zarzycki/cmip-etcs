#!/bin/bash

THISDIR=$PWD
FZRADIR="/glade/derecho/scratch/zarzycki/FZRA/"

# Define years to process
START_YEAR=1980
END_YEAR=2020

cd $FZRADIR/ptypes/ || { echo "Failed to change directory to $FZRADIR/ptypes/. Exiting."; exit 1; }
pwd

# Define datasets and their input file patterns
declare -A datasets=(
    ["JRA"]="JRA/JRA_PTYPE_YYYY.nc"
    ["CFSR"]="CFSR/CFSR_PTYPE_YYYY.nc"
    ["ECMWF"]="ECMWF_internal/ecmwf_reformat_YYYY.nc"
    ["ERA5"]="ERA5/ERAML_YYYY.nc"
    ["CR20"]="CR20/CR20_PTYPE_YYYY.nc"
)

# Loop through each dataset
for dataset in "${!datasets[@]}"; do
    echo "Processing ${dataset} dataset..."

    # Get the file pattern for this dataset
    pattern="${datasets[$dataset]}"

    # Loop through years
    for year in $(seq $START_YEAR $END_YEAR); do
        # Replace YYYY with actual year in file pattern
        input_file="${pattern/YYYY/$year}"
        output_file="${dataset}-${year}.nc"

        echo "Processing ${dataset} data for year ${year}..."

        # Check if input file exists
        if [ -f "$input_file" ]; then
            if [ "$dataset" = "ECMWF" ]; then
                # Special processing for ECMWF
                ncremap -v ptype -a neareststod -d $THISDIR/target.nc -i "$input_file" -o "$output_file"
                if [ $? -eq 0 ]; then
                    # Additional ECMWF processing
                    ncks -O -d time,,,6 "$output_file" "$output_file"
                    ncrename -O -v ptype,PTYPE "$output_file" "$output_file"
                    echo "Successfully processed ${dataset} ${year} with ECMWF-specific steps"
                else
                    echo "Error processing ${dataset} ${year}"
                fi
            elif [ "$dataset" = "ERA5" ]; then
                # Special processing for ERA5
                ncremap -v PTYPE -a neareststod -d $THISDIR/target.nc -i "$input_file" -o "$output_file"
                if [ $? -eq 0 ]; then
                    # Additional ERA5 processing
                    ncks -O -d time,,,6 "$output_file" "$output_file"
                    echo "Successfully processed ${dataset} ${year} with ERA5-specific steps"
                else
                    echo "Error processing ${dataset} ${year}"
                fi
            else
                # Standard processing for other datasets
                ncremap -v PTYPE -a neareststod -d $THISDIR/target.nc -i "$input_file" -o "$output_file"
                if [ $? -eq 0 ]; then
                    echo "Successfully processed ${dataset} ${year}"
                else
                    echo "Error processing ${dataset} ${year}"
                fi
            fi
        else
            echo "Warning: Input file not found: ${input_file}"
        fi
    done
done
