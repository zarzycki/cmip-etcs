#!/bin/bash -l
#
# ==============================================================================
#  ETC TRACKING WITH TEMPESTEXTREMES -- ANNOTATED STARTER SCRIPT
# ==============================================================================
#
#  WHAT THIS DOES
#  --------------
#  Finds extratropical cyclones (ETCs) in ERA5 reanalysis and stitches them into
#  storm tracks, restricted to storms that pass over/near the continental US.
#
#  It is deliberately small: a few months of data, so it finishes in minutes on
#  a login node. The production scripts in the parent directory
#  (../etc-tracking-cmip6.sh, ../etc-tracking-era5.sh) do the same thing over
#  40+ years and many models; read them AFTER you understand this one.
#
#  THE TWO-STEP IDEA
#  -----------------
#  Nearly all TempestExtremes feature tracking is two stages:
#
#    1. DetectNodes  -- Looks at each time slice INDEPENDENTLY and asks
#                       "where are the candidate features right now?"
#                       Output: a text file of candidate points, no time
#                       continuity, no notion of a "storm".
#
#    2. StitchNodes  -- Reads that candidate file and connects points ACROSS
#                       TIME into trajectories, applying rules about how far a
#                       storm can move, how long it must last, etc.
#                       Output: a text file of storm tracks.
#
#  Everything downstream (compositing, precip attribution, statistics) operates
#  on the stage-2 trajectory file. See 02-nodefile-format.md for what these
#  text files actually look like -- read that doc, it is the main thing that
#  makes TempestExtremes click.
#
#  HOW TO RUN
#  ----------
#    conda activate tempest-extremes
#    ./01-detect-stitch-era5.sh
#
# ==============================================================================


# Stop immediately if any command fails, or if an undefined variable is used.
# Without this, a failed DetectNodes would silently produce an empty file and
# StitchNodes would "succeed" with zero storms -- a very confusing way to lose
# an afternoon.
set -euo pipefail


# ============================== USER OPTIONS ==================================

# Which conda environment holds the TempestExtremes binaries. If `which
# DetectNodes` already works you can ignore this.
TEMPESTDIR=/glade/work/zarzycki/conda-envs/tempest-extremes/bin

# Where ERA5 lives on the NCAR campaign store (RDA dataset ds633.0).
# Mean sea level pressure, 0.25 degree, HOURLY, one file per month.
ERA5DIR=/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc

# Time period. Start with a single winter season -- roughly 45 s of compute per
# month. Once this works, try MONTHS="01 02 03 04 05 06 07 08 09 10 11 12" for a
# full year (~10 min), and then look at ../etc-tracking-era5.sh for how the
# multi-decade production runs are structured (hint: they need a PBS job).
YEAR=1990
MONTHS="01 02 03"

# Where to put results.
OUTDIR="./output"

# =============================== END OPTIONS ==================================


# Put the TempestExtremes binaries on PATH for this script only.
export PATH="${TEMPESTDIR}:${PATH}"

mkdir -p "${OUTDIR}"

# Names for the two stage outputs.
CANDIDATEFILE="${OUTDIR}/candidates_ERA5_${YEAR}.txt"
TRAJFILE="${OUTDIR}/tracks_ERA5_${YEAR}.txt"
FILELIST="${OUTDIR}/filelist_ERA5_${YEAR}.txt"

# Scratch prefix used while DetectNodes runs (explained in step 1).
CANDPREFIX="${OUTDIR}/candpart_${YEAR}"


# ==============================================================================
#  STEP 0: BUILD THE INPUT FILE LIST
# ==============================================================================
#
#  TempestExtremes can take a single file (--in_data) or a text file listing
#  many files, one per line (--in_data_list). The list form is what you want
#  for anything real, because the tracker needs to follow storms across the
#  boundaries between monthly files.
#
#  GOTCHA: use `find`, not `ls`. On most NCAR accounts `ls` is aliased to add
#  color, which injects invisible ANSI escape codes into the file. TempestExtremes
#  then reports "Cannot open input file" for a path that looks perfectly correct
#  on your screen. This bites everyone once.

echo "=== Building file list ==="
rm -f "${FILELIST}"
for MM in ${MONTHS}; do
  find "${ERA5DIR}/${YEAR}${MM}/" -name "e5.oper.an.sfc.128_151_msl.ll025sc.*.nc"
done | sort > "${FILELIST}"

echo "Found $(wc -l < "${FILELIST}") input file(s):"
cat "${FILELIST}"


# ==============================================================================
#  STEP 1: DETECTNODES -- FIND CANDIDATE CYCLONE CENTERS
# ==============================================================================
#
#  An ETC is identified here as a local minimum in sea level pressure that sits
#  inside a CLOSED PRESSURE CONTOUR. The closed-contour test is what separates a
#  real cyclone from a trivial dimple in the pressure field: we demand that
#  pressure rise by a set amount in every direction as you walk away from the
#  minimum.
#
#  The parameter choices below follow the configuration used in this project's
#  published work (Zarzycki 2018; Gore et al. 2023). They are NOT universal
#  constants -- they are a tuned definition of "what counts as an ETC", and part
#  of your job on this project will be defending or revisiting them.

MSLNAME="MSL"        # Variable name in the ERA5 netCDF files. In CMIP6 it is
                     # "psl"; in CESM output it is "PSL". Same physical field,
                     # three different names -- a recurring annoyance.

CONTOUR_DELTA=200.0  # Pa. Pressure must INCREASE by 200 Pa (2 hPa) from the
                     # minimum outward for the contour to count as closed.
                     # Larger => only stronger, better-defined lows survive.

CONTOUR_DIST=6.0     # Great-circle degrees. The 200 Pa rise must be achieved
                     # within 6 degrees of the center. This sets the SPATIAL
                     # SCALE of the feature: 6 deg is synoptic (ETC-sized).
                     # TC tracking uses much smaller values.

MERGEDIST=6.0        # Great-circle degrees. If two candidate minima are within
                     # 6 degrees of each other, keep only the stronger one.
                     # Prevents a single broad cyclone with a ragged pressure
                     # field from being reported as three separate storms.

TIMESTRIDE=6         # ERA5 is hourly. Take every 6th time slice => 6-hourly.
                     # This matches the CMIP6 output frequency, so ERA5 and the
                     # models are compared on equal footing. See
                     # ../tracker-dt-sens/ for the experiments on how much this
                     # choice matters -- it matters more than you would think.

# Detection domain. Restricting the search box makes this tutorial fast.
#
# IMPORTANT SUBTLETY: this is NOT how you should subset for real science. The
# production runs detect GLOBALLY and then filter trajectories afterwards (see
# the --threshold in step 2). Detecting inside a box truncates any storm that
# crosses the box edge, so its track appears to begin or end at the boundary.
# The box here is drawn generously around the North Atlantic/North America so
# that storms reaching the US still have most of their life cycle inside it.
MINLAT=15.0
MAXLAT=75.0
MINLON=180.0
MAXLON=350.0

echo ""
echo "=== STEP 1: DetectNodes ==="


# A wrinkle you need to know about: when given --in_data_list, DetectNodes
# writes ONE OUTPUT FILE PER INPUT FILE, named <out>000000.dat, <out>000001.dat,
# and so on -- it does NOT produce the single file you named in --out. So we
# point --out at a temporary prefix and concatenate the pieces ourselves below.
#
# Every production script in this repo does the same dance. If you ever see
# "Unable to open input file" from StitchNodes, this is why.
#
# (There is also an --out_file_list option that records the piece names, but it
# cannot be combined with --out, so globbing is the simpler route.)

rm -f "${CANDPREFIX}"*.dat

DetectNodes \
  --in_data_list "${FILELIST}" \
  --out "${CANDPREFIX}" \
  --searchbymin "${MSLNAME}" \
  --closedcontourcmd "${MSLNAME},${CONTOUR_DELTA},${CONTOUR_DIST},0" \
  --mergedist "${MERGEDIST}" \
  --outputcmd "${MSLNAME},min,0" \
  --timestride "${TIMESTRIDE}" \
  --latname "latitude" \
  --lonname "longitude" \
  --minlat "${MINLAT}" --maxlat "${MAXLAT}" \
  --minlon "${MINLON}" --maxlon "${MAXLON}" \
  --verbosity 0

# A note on the two commands above that look similar but are not:
#
#   --closedcontourcmd MSL,200.0,6.0,0
#         A *test*. Points failing it are thrown away.
#         Format: var,delta,dist,minmaxdist
#
#   --outputcmd MSL,min,0
#         A *record*. For every point that survived, write the minimum MSL
#         within 0 degrees (i.e. the value at the point itself) as an extra
#         column. This is how you attach storm properties to a track.
#         Format: var,operation,dist
#
# Want the wind speed at each cyclone center too? Add another --outputcmd entry
# (separated by semicolons) and add a matching name to --in_fmt in step 2. That
# is exactly the mechanism you will use for the wind metrics this project needs.

# IMPORTANT: TempestExtremes prints its errors as "EXCEPTION ..." but still
# exits with status 0, so `set -e` will NOT catch a failed run. Check for real
# output explicitly, or you will spend a long time debugging an empty plot.
if ! compgen -G "${CANDPREFIX}*.dat" > /dev/null; then
  echo "ERROR: DetectNodes produced no output. Scroll up for an EXCEPTION line." >&2
  exit 1
fi

# Glue the per-file pieces into one candidate file, in time order. The numeric
# suffix is zero-padded, so a plain sort puts them in the right order.
rm -f "${CANDIDATEFILE}"
for PART in $(ls -1 "${CANDPREFIX}"*.dat | sort); do
  cat "${PART}" >> "${CANDIDATEFILE}"
done

# Tidy up the intermediate pieces.
rm -f "${CANDPREFIX}"*.dat

echo "DetectNodes wrote $(wc -l < "${CANDIDATEFILE}") lines to ${CANDIDATEFILE}"


# ==============================================================================
#  STEP 2: STITCHNODES -- CONNECT CANDIDATES INTO TRACKS
# ==============================================================================
#
#  Now we impose time continuity. StitchNodes walks forward in time and links
#  candidate points into paths, then throws away paths that do not behave like
#  real ETCs.

RANGE=6.0            # Great-circle degrees. Maximum distance a storm may travel
                     # between consecutive (6-hourly) time steps. 6 deg / 6 h is
                     # about 110 km/h of translation -- generous but physical.

MINTIME="60h"        # A track must last at least 60 hours (2.5 days) to be kept.
                     # Removes transient pressure wiggles that briefly satisfy
                     # the closed-contour test.
                     #
                     # NOTE: older scripts in this repo pass `--minlength 10`
                     # (10 time steps). That flag NO LONGER EXISTS in the current
                     # conda-forge TempestExtremes build and will cause an error.
                     # Use the duration form. This is one reason the old
                     # ../etc-tracking-cmip6.sh will not run as-is today.

MAXGAP="18h"         # A track may skip up to 18 h (3 steps) without a detected
                     # center and still be considered the same storm. Cyclones
                     # can temporarily lose their closed contour while remaining
                     # coherent, e.g. during occlusion or over steep terrain.

MINENDPOINTDIST=12.0 # Great-circle degrees between first and last point. This
                     # discards quasi-stationary features -- persistent thermal
                     # or orographic lows that sit in one place for days and are
                     # not travelling storms.

# Geographic filter, expressed on the trajectory rather than on detection.
#
# Read as: latitude > 24 AND longitude > 234 AND latitude < 52 AND longitude
# < 294. Longitudes are 0-360, so this is roughly 126W-66W, 24N-52N: CONUS.
#
# The trailing ",1" is the count, and it is the important part: the condition
# must hold for AT LEAST 1 point in the track. So a storm that spins up over the
# Pacific, crosses the US, and dies over the Atlantic is kept IN ITS ENTIRETY,
# including the portions outside the box. You are selecting whole storms that
# visit the region, not clipping tracks to the region.
#
# Raise the count (e.g. ",4") to require a storm to linger over the domain.
CONUSFILTER="lat,>,24,1;lon,>,234,1;lat,<,52,1;lon,<,294,1"

# --in_fmt tells StitchNodes what the columns of the candidate file MEAN.
# It must list the --outputcmd variables from step 1, in order. We recorded one
# extra variable (the MSL minimum), which we name "slp" here.
#
# If you add outputs in step 1, you MUST extend this string or every downstream
# tool will silently read the wrong column.
INFMT="lon,lat,slp"

echo ""
echo "=== STEP 2: StitchNodes ==="

StitchNodes \
  --in "${CANDIDATEFILE}" \
  --out "${TRAJFILE}" \
  --in_fmt "${INFMT}" \
  --range "${RANGE}" \
  --mintime "${MINTIME}" \
  --maxgap "${MAXGAP}" \
  --min_endpoint_dist "${MINENDPOINTDIST}" \
  --threshold "${CONUSFILTER}"

if [ ! -s "${TRAJFILE}" ]; then
  echo "ERROR: StitchNodes produced no tracks. Scroll up for an EXCEPTION line." >&2
  exit 1
fi

NTRACKS=$(grep -c '^start' "${TRAJFILE}" || true)

echo ""
echo "=============================================="
echo " Done."
echo " Candidates : ${CANDIDATEFILE}"
echo " Tracks     : ${TRAJFILE}"
echo " ETCs found : ${NTRACKS}"
echo "=============================================="
echo ""
echo "Next: read 02-nodefile-format.md, then run"
echo "  python 03-plot-tracks.py ${TRAJFILE}"
