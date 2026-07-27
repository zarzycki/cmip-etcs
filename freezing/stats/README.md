# Stats and verification scripts

These scripts compare derived precip / precip-type products against a reference dataset. Both are configured by editing the constants in their `if __name__ == "__main__":` block (no CLI args), and both are typically run in the background with output captured to a log file, e.g.:

```
nohup python metrics-precip-corr-rmse.py > precip.out &
nohup python metrics-ptype-flags.py > program.out &
```

`*.out` logs are gitignored — they're just run records for a particular scratch dataset/path and aren't meant to be committed.

### metrics-precip-corr-rmse.py

Verifies gridded precipitation against a reference product.

- Reference dataset: `IMERG`
- Test datasets: `ERA5`, `JRA`, `CFSR`, `CR20`, `TEST`
- Years: 2001-2015
- Time offsets tested: `-6, 0, 6` hours (to check sensitivity to timestamp alignment)
- Data directory: hardcoded `FZRAPATH` (currently `/glade/derecho/scratch/zarzycki/FZRA/precip/`)

For each dataset/year/offset it computes RMSE, Pearson correlation, bias, POD, FAR, and CSI (the latter three using a 0.1 mm/day threshold), then prints a final summary table across all years.

### metrics-ptype-flags.py

Verifies gridded precipitation-*type* classification (rain/freezing rain/ice pellets/snow, encoded as class labels `1, 2, 4, 8`) against a reference product.

- Reference dataset: `ECMWF`
- Comparison datasets: `ERA5`, `CFSR`, `CR20`, `JRA`
- Years: 1980-2016
- Data directory: hardcoded `FZRAPATH` (currently `/glade/derecho/scratch/zarzycki/FZRA/ptypes/`)

For each dataset/year it computes overall accuracy plus per-class precision/recall/F1, then prints a mean ± std summary across all years per dataset.
