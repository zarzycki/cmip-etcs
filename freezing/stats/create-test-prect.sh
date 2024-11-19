#!/bin/bash

FZRADIR="/glade/derecho/scratch/zarzycki/FZRA/"

cd $FZRADIR/precip/ || { echo "Failed to change directory to $FZRADIR/precip/. Exiting."; exit 1; }
pwd

# Remove any existing symlinked files
rm -fv pre-TEST-*.nc

# Function to check if a year is leap year
is_leap_year() {
    year=$1
    if (( year % 4 == 0 )) && (( year % 100 != 0 )) || (( year % 400 == 0 )); then
        return 0  # true
    else
        return 1  # false
    fi
}

# Create arrays for leap years and non-leap years
declare -a leap_years=()
declare -a non_leap_years=()

# Sort years into appropriate arrays
for file in pre-ERA5-*.nc; do
    year=$(echo $file | grep -o '[0-9]\{4\}')
    if is_leap_year $year; then
        leap_years+=($year)
    else
        non_leap_years+=($year)
    fi
done

# Sort arrays numerically
IFS=$'\n' leap_years=($(sort -n <<<"${leap_years[*]}"))
IFS=$'\n' non_leap_years=($(sort -n <<<"${non_leap_years[*]}"))

# Create symbolic links for non-leap years (+1 offset)
for ((i=0; i<${#non_leap_years[@]}; i++)); do
    source_year=${non_leap_years[$i]}
    # Get target year (wrap around if needed)
    next_idx=$(( (i + 1) % ${#non_leap_years[@]} ))
    target_year=${non_leap_years[$next_idx]}
    ln -s pre-ERA5-${source_year}.nc pre-TEST-${target_year}.nc
done

# Create symbolic links for leap years (+4 offset)
for ((i=0; i<${#leap_years[@]}; i++)); do
    source_year=${leap_years[$i]}
    # Get target year (wrap around if needed)
    next_idx=$(( (i + 1) % ${#leap_years[@]} ))
    target_year=${leap_years[$next_idx]}
    ln -s pre-ERA5-${source_year}.nc pre-TEST-${target_year}.nc
done
