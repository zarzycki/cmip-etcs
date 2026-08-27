# Getting started: tracking extratropical cyclones with TempestExtremes

A self-contained introduction to the ETC tracking pipeline used in this project.
Assumes you are comfortable with gridded netCDF data but have never used
TempestExtremes. Budget an afternoon.

Everything here runs on Derecho in a few minutes. The production scripts in the
parent directory do the same thing over 40+ years and dozens of models — read
those *after* you have worked through this.

---

## Project context

This work is part of an ESTCP-funded project, *Informing DoD Infrastructure
Investments Through Storm-Typing Climatologies and Projections* (see
`../NARRATIVE.txt`). The idea: most extreme US precipitation and wind comes from
four storm types — tropical cyclones, **extratropical cyclones**, mesoscale
convective systems, and atmospheric rivers. If you can detect and track each
type in climate model output, you can attribute hazards to storm type, judge
which models get them right, and project how they change.

There are four technical objectives:

| | | |
|---|---|---|
| **TO1** | Track all four storm types in CMIP6 and high-resolution regional models | ← you are here |
| **TO2** | Evaluate against observations and reanalysis; score models by region | |
| **TO3** | Project changes at 2035, 2065, 2100 under SSP2-4.5 and SSP5-8.5 | |
| **TO4** | Build the Coastal Storm Hazard Mapper decision-support tool | (other institutions) |

**Our group leads the ETC piece of TO1, TO2, and TO3.** So the through-line for
your work is: track ETCs reliably → attribute precipitation and wind to them →
show how well models reproduce that → show how it changes in the future.

This tutorial is the first step of the first objective.

---

## Setup

```bash
conda activate tempest-extremes
```

That environment has the TempestExtremes binaries plus python, xarray, numpy,
matplotlib, cartopy, cdo, and the NCO tools. Check it works:

```bash
which DetectNodes StitchNodes
```

Both should resolve. If not, the binaries are in
`/glade/work/zarzycki/conda-envs/tempest-extremes/bin`.

---

## Work through these in order

### 1. Run the tracker — `01-detect-stitch-era5.sh`

```bash
./01-detect-stitch-era5.sh
```

Roughly 2.5 minutes for three months of ERA5. Finds ETCs and stitches them into
tracks over the continental US.

**Read the script, don't just run it.** It is commented far more heavily than
anything you would normally write, and the comments are the actual content — what
each tracking parameter means physically and why it has the value it has.

### 2. Understand the output — `02-nodefile-format.md`

The single most important document here. TempestExtremes passes data around as
ASCII text files with positionally-defined columns, which is unfamiliar if you
are coming from gridded netCDF. Nearly every confusing error downstream traces
back to misunderstanding these files.

### 3. Look at the storms — `03-plot-tracks.py`

```bash
python 03-plot-tracks.py output/tracks_ERA5_1990.txt --out output/etc_tracks.png
```

Plots the tracks and a track-density map. Contains a trajectory-file parser you
should copy into your own analysis code — there is no library for this, everyone
writes it once.

### 4. Verify and experiment — `04-check-your-work.md`

Compare against the reference output in `reference/`, then deliberately break
the tracking parameters one at a time to build intuition for what each controls.
This part matters more than it sounds like it does.

---

## The mental model

Two stages, always:

```
   ERA5 / CMIP6 netCDF
            |
            |  DetectNodes    "where are candidate cyclones at each time?"
            v
   candidate file (ASCII)     independent points, no time continuity
            |
            |  StitchNodes    "which points belong to the same storm?"
            v
   trajectory file (ASCII)    storm tracks
            |
            +--> HistogramNodes   track density
            +--> NodeFileFilter   attribute precip/wind to storms   <- the science
            +--> NodeFileEditor   subset tracks by intensity
            +--> NodeFileCompose  storm-centred composites
```

`NodeFileFilter` is where storm-typing actually happens: mask a precipitation
field to only the precipitation falling near tracked ETCs, and you have
ETC-attributed precipitation. Repeat for each storm type and you can decompose
total rainfall by what produced it.

---

## Scaling up

The conda environment ships the **serial** TempestExtremes build
(`tempest-extremes-2.4.2-nompi_*`). Don't launch it under `mpiexec` — every rank
would process the whole file list and clobber the same output prefix, silently.

That is fine, because `DetectNodes` is embarrassingly parallel: its MPI mode just
round-robins entries of `--in_data_list` across ranks, with no communication. To
go faster, split the file list into N chunks, run one serial `DetectNodes` per
chunk as a PBS job array, and `cat` the candidate files together — which is the
same concatenation step `01-detect-stitch-era5.sh` already does. `StitchNodes`
then runs once over the combined file.

MPI builds do exist upstream
(`conda install -c conda-forge 'tempest-extremes=*=mpi_mpich_*'`), but they use
conda MPICH rather than Cray MPICH, so they won't use Derecho's interconnect,
and they pull an MPI netCDF/HDF5 stack that can disturb the python side of the
environment. Not worth it for this workload.

Note also that `../etc-tracking-era5.sh` calls `mpiexec_mpt`, the Cheyenne
SGI-MPT launcher, which does not exist on Derecho.

---

## What's in the parent directory

Orientation, roughly in order of usefulness to you:

| Path | What it is |
|---|---|
| `../etc-tracking-era5.sh` | Production ERA5 tracking, full record |
| `../etc-tracking-cmip6.sh` | Production CMIP6 tracking, per-model config blocks |
| `../etc-tracking-lens.sh` | CESM1 Large Ensemble variant |
| `../trajs/` | Already-tracked output: 13 models + 7 reanalyses |
| `../tracker-dt-sens/` | Experiments on hourly vs 6-hourly tracking sensitivity |
| `../data-process/` | Topography file generation, time-coordinate fixes |
| `../freezing/` | Precipitation-type (rain/snow/freezing rain) algorithms |
| `../for-gmd-v2/` | Figures and analysis for a methods paper |
| `../plotting/` | Older NCL plotting scripts |

Two caveats on the older material:

- **The production scripts predate Derecho.** They contain hostname branches for
  Cheyenne, ACI, and a laptop, and reference `/glade/scratch`, which no longer
  exists. Paths need updating before they will run.
- **`--minlength` was removed** from TempestExtremes. Scripts using it will
  error against the current build; use `--mintime` instead.

Do not try to fix all of that up front. Get this tutorial working first.

---

## Reference

- TempestExtremes docs: https://climate.ucdavis.edu/tempestextremes.php
- Ullrich & Zarzycki (2017), *TempestExtremes v1.0*, Geosci. Model Dev.
- Ullrich et al. (2021), *TempestExtremes v2.1*, Geosci. Model Dev.
- Zarzycki (2018) — ETC tracking applied to Regional Snowfall Index
- Reed et al. (2023) — the storm-typing precipitation decomposition this project builds on

Run any TempestExtremes tool with no arguments to print its full argument list.
That is often faster than the documentation.
