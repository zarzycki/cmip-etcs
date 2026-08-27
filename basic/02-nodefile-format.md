# The node file format (read this one)

TempestExtremes passes data between its tools as **plain ASCII text files**, not
netCDF. If you are used to gridded data this is the unfamiliar part, and almost
every confusing error message downstream traces back to a misunderstanding of
these files. They are worth ten minutes.

There are two of them, and they look similar but mean different things.

---

## 1. The candidate file (output of `DetectNodes`)

Independent snapshots in time. No storms yet — just "here are the points that
passed the test at this instant."

```
1990	1	1	7	0
	864	84	216.000000	69.000000	1.016409e+05
	1269	118	317.250000	60.500000	9.727562e+04
	919	122	229.750000	59.500000	9.895738e+04
1990	1	1	7	6
	871	85	217.750000	68.750000	1.015883e+05
	...
```

**Header lines** (no leading tab) start each time slice:

```
1990     1       1        7         0
YEAR    MONTH   DAY   n_points    HOUR
```

Careful — the 4th field is the **number of candidate points at this time**, not
part of the date. The hour is last.

**Data lines** (leading tab), one per candidate:

```
	864	84	216.000000	69.000000	1.016409e+05
	 i	 j	     lon	     lat	      MSL
```

- `i`, `j` — the grid indices the point was found at. For a rectilinear grid
  like ERA5 these are the longitude and latitude indices. For unstructured grids
  there is only one index and `j` is absent. You rarely need these.
- `lon`, `lat` — the position in degrees. **Longitude is 0–360**, not −180–180.
  This trips people up constantly when plotting or comparing against other
  datasets.
- Everything after that comes from your `--outputcmd` entries, **in the order
  you listed them**. Here we asked for one thing (`MSL,min,0`) so there is one
  extra column: 1.016409e+05 Pa = 1016.4 hPa.

If you had asked for `--outputcmd "MSL,min,0;VAR_10U,max,2.0"` there would be two
extra columns and the second would be max 10 m wind within 2 degrees.

---

## 2. The trajectory file (output of `StitchNodes`)

Same points, now grouped into storms with time continuity.

```
start	29	1990	1	1	0
	1131	174	282.750000	46.500000	9.938538e+04	1990	1	1	0
	1144	169	286.000000	47.750000	9.895800e+04	1990	1	1	6
	1160	165	290.000000	48.750000	9.846319e+04	1990	1	1	12
	1189	161	297.250000	49.750000	9.749319e+04	1990	1	1	18
start	17	1990	1	4	6
	...
```

**Header lines** begin with the literal word `start`:

```
start     29     1990    1      1     0
        npts     YEAR   MON    DAY   HOUR     <- of the FIRST point in the track
```

`npts` is how many data lines follow before the next `start`.

**Data lines** are the candidate-file columns **plus a full date stamp appended
to each point**:

```
	1131	174	282.750000	46.500000	9.938538e+04	1990	1	1	0
	   i	  j	     lon	     lat	         slp	YEAR	MON	DAY	HR
```

That per-point timestamp is the whole difference. Now you can ask when a storm
was deepest, how fast it moved, where it was at landfall.

So the example above is a cyclone at 46.5N, 282.75E (= 77.25W, near Lake
Ontario) at 00Z on 1 Jan 1990 with a central pressure of 993.9 hPa, deepening to
974.9 hPa by 18Z as it tracks northeast. That is a textbook East Coast winter
storm.

---

## 3. `--in_fmt`: the part that actually matters

Nothing in the file says what the columns are. Every downstream tool is told
separately, via `--in_fmt`:

```
--in_fmt "lon,lat,slp"
```

This names the columns **after** `i`/`j` and **before** the date stamp. The
names are yours to choose (`slp` is a convention, not a keyword) — but they must
match, in order, the `--outputcmd` list you gave `DetectNodes`.

Three consequences worth internalising:

1. **`--in_fmt` is not validated.** Get the order wrong and tools will happily
   read pressure as latitude. There is no error, just nonsense results.
2. **Every tool in the chain needs it.** `StitchNodes`, `NodeFileFilter`,
   `NodeFileEditor`, `NodeFileCompose` all take `--in_fmt`, and it must be
   consistent across all of them.
3. **Adding a variable means editing several places at once**: the
   `--outputcmd` in `DetectNodes`, and the `--in_fmt` in every subsequent call.

You can see this in the production scripts: `../etc-tracking-cmip6.sh` uses
`--in_fmt "lon,lat,slp,phis"` because it records surface geopotential
(`PHIS,max,0`) alongside pressure, so it can tell whether a "cyclone" is really
a terrain artifact over the Rockies. Our tutorial script has no `phis` column,
which is why its `--in_fmt` is shorter.

---

## 4. Where this goes next

Once you have a trajectory file, the rest of the TempestExtremes toolkit reads
it and does the science:

| Tool | What it does | Used in this project for |
|---|---|---|
| `HistogramNodes` | Bins track points onto a lat/lon grid | Track density maps, model-vs-ERA5 comparison |
| `NodeFileFilter` | Masks a gridded field to points near tracks | Attributing precipitation/wind to ETCs (TO1) |
| `NodeFileEditor` | Filters/edits tracks by column value | Pulling out the `slp <= 990 hPa` strong subset |
| `NodeFileCompose` | Storm-centered composites of a gridded field | Mean ETC precip/wind structure |

`NodeFileFilter` is the one that connects storm tracks back to hazards, and it
is the heart of the storm-typing method this whole project is built on: take a
precipitation field, keep only the precipitation that fell within some radius of
a tracked ETC, and you have "ETC-attributed precipitation." Do that for TCs,
ETCs, MCSs, and ARs and you can decompose total rainfall by storm type — that is
what Figure 1 of the proposal shows.

---

## 5. Reading these files in Python

There is no library for this; everyone writes a small parser. `03-plot-tracks.py`
in this directory has one you can copy — the logic is just "lines starting with
`start` (or lacking a leading tab) open a new storm; everything else is a point
in the current storm."
