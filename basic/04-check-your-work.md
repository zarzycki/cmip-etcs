# Did it work?

## The expected answer

Running `01-detect-stitch-era5.sh` unmodified (ERA5, Jan–Mar 1990) should give
you exactly this:

```
Storms            : 19
Mean duration     : 115.6 h (min 60, max 294)
Mean min pressure : 978.1 hPa (deepest 932.9)
Total storm points: 349
```

and a trajectory file of 368 lines with MD5 `586b3d133c974aa824f65e5457346d37`.

A byte-identical reference copy lives in `reference/era5_1990JFM_tracks.txt.gz`.
Check yourself against it directly:

```bash
zcat reference/era5_1990JFM_tracks.txt.gz | diff - output/tracks_ERA5_1990.txt \
  && echo "MATCH"
```

TempestExtremes is deterministic, so this should match exactly. If it does not,
something in your environment or input data differs — worth chasing down before
you build anything on top of it.

## Does the answer make physical sense?

Reproducing a checksum only proves you ran the same code. Confirm the result is
*meteorologically* believable, which is the habit that actually matters:

- **19 storms in three winter months over CONUS** — plausible. Winter is the
  active ETC season, and the filters (60 h minimum lifetime, 12° minimum
  displacement) are deliberately strict, so this is a subset of "every low that
  ever appeared", not a full count.
- **Mean minimum pressure 978 hPa** — reasonable for storms that survive those
  filters. The 932.9 hPa outlier is a genuine bomb cyclone; find it and look it
  up.
- **In the tracks plot**, look for the classic North American storm paths:
  Alberta Clippers dropping southeast out of Canada, Colorado lows forming in
  the Rockies' lee, and East Coast storms intensifying as they run up toward the
  Canadian Maritimes. If your tracks do not look like that, something is wrong
  regardless of what the checksum says.
- **In the density plot**, the maxima should sit downstream of the Rockies and
  off the eastern seaboard, not scattered randomly.

## Things worth breaking on purpose

The parameters are the science, so develop intuition for what each one does by
changing exactly one at a time and re-running:

| Change | In file | What you should see |
|---|---|---|
| `CONTOUR_DELTA` 200 → 400 | `01-...sh` | Far fewer storms; only deep, well-defined lows pass |
| `MINTIME` `60h` → `24h` | `01-...sh` | Many more, shorter-lived tracks |
| `MINENDPOINTDIST` 12.0 → 0.0 | `01-...sh` | Quasi-stationary lows reappear, especially over terrain |
| `MERGEDIST` 6.0 → 1.0 | `01-...sh` | One storm starts fragmenting into several |
| `TIMESTRIDE` 6 → 1 | `01-...sh` | Hourly tracking — slower, and the track count changes. See `../tracker-dt-sens/` |
| CONUS `--threshold` count `,1` → `,8` | `01-...sh` | Only storms that linger ≥ 2 days over the domain |

Keep notes on what you find. Justifying these choices is a real part of the
project, not busywork — the numbers currently in the production scripts came
from exactly this kind of exercise.

## Common failure modes

**"Cannot open input file" but the path looks right.**
Colour codes from an aliased `ls` are in your file list. Use `find`.

**"Unable to open input file <candidates>" from StitchNodes.**
With `--in_data_list`, DetectNodes writes `<out>000000.dat`, `<out>000001.dat`,
… rather than the single file you named. Concatenate them first (the script
does this in step 1).

**"Invalid argument --minlength".**
Old flag, removed from the current TempestExtremes build. Use `--mintime` with a
duration like `60h`. Several of the older scripts in the parent directory still
use it and will need updating.

**Zero storms found, no error message.**
TempestExtremes prints `EXCEPTION` and still exits 0, so `set -e` and `&&`
chains will not catch failures. Always scroll up and look. This is also why the
script explicitly tests for non-empty output.

**Everything runs but the plot is blank/nonsense.**
Your `--in_fmt` probably does not match the `--outputcmd` list. Nothing
validates this; columns are read positionally.

## When this all makes sense

Move on to a CMIP6 model rather than another reanalysis. That is where the real
complications live, and it is the actual TO1 deliverable:

- The pressure variable is `psl`, not `MSL`, and dimensions are `lat`/`lon`, not
  `latitude`/`longitude`.
- You need a topography file so you can record `PHIS` alongside pressure and
  screen out spurious terrain-driven "cyclones" — see
  `../data-process/make_topo_file.ncl`.
- Several models have broken or non-standard time coordinates (the notes in
  `../etc-tracking-cmip6.sh` flag NESM3 and BCC-CSM2-MR).
- Model calendars are often `noleap` or `360_day`, which StitchNodes needs told
  via `--caltype`.

Then the science question this project actually asks: does the model produce the
same ETC frequency, intensity distribution, and track density as ERA5? That
comparison is TO2.
