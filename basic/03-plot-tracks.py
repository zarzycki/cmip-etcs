#!/usr/bin/env python
"""
Plot ETC tracks and track density from a TempestExtremes trajectory file.

Usage:
    conda activate tempest-extremes
    python 03-plot-tracks.py output/tracks_ERA5_1990.txt
    python 03-plot-tracks.py output/tracks_ERA5_1990.txt --out mytracks.png

The parser at the top (`read_trajectories`) is the part worth stealing for your
own analysis scripts. There is no standard library for reading these files, so
everyone ends up writing this once. See 02-nodefile-format.md for what the
format means.
"""

import argparse
import gzip
from datetime import datetime

import numpy as np
import matplotlib

matplotlib.use("Agg")  # Render to file; no display on Derecho login nodes.
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# ---------------------------------------------------------------------------
# Reading the trajectory file
# ---------------------------------------------------------------------------

def read_trajectories(filename, in_fmt="lon,lat,slp"):
    """Parse a TempestExtremes StitchNodes trajectory file.

    The file looks like:

        start   29      1990    1       1       0
                1131    174     282.75  46.50   9.938538e+04    1990  1  1  0
                1144    169     286.00  47.75   9.895800e+04    1990  1  1  6
        start   17      1990    1       4       6
                ...

    Header lines begin with the word 'start'; everything else is a point
    belonging to the most recently opened track.

    Each data line is:  i, j, <in_fmt columns...>, year, month, day, hour

    Parameters
    ----------
    filename : str
        Path to the trajectory file. Plain text or gzipped (.gz) both work.
    in_fmt : str
        The same string passed to StitchNodes --in_fmt. Must list the columns
        that sit between the (i, j) grid indices and the trailing date stamp.

    Returns
    -------
    list of dict
        One dict per storm, with a numpy array for each named column plus a
        'time' array of datetime objects.
    """
    columns = in_fmt.split(",")

    # Layout: 2 index columns, then the named columns, then 4 date columns.
    n_leading = 2
    n_date = 4

    opener = gzip.open if str(filename).endswith(".gz") else open

    tracks = []
    current = None

    with opener(filename, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue

            if line.lstrip().startswith("start"):
                # New storm begins. We deliberately ignore the header's own
                # date/count fields and rebuild them from the data lines, which
                # is more robust than trusting the header.
                current = {name: [] for name in columns}
                current["time"] = []
                tracks.append(current)
                continue

            if current is None:
                # Data before any 'start' line -- this is probably a DetectNodes
                # candidate file, not a StitchNodes trajectory file.
                raise ValueError(
                    f"{filename}: found a data line before any 'start' header. "
                    "Is this a candidate file rather than a trajectory file?"
                )

            parts = line.split()
            expected = n_leading + len(columns) + n_date
            if len(parts) != expected:
                raise ValueError(
                    f"{filename}: expected {expected} fields for "
                    f"in_fmt='{in_fmt}' but found {len(parts)}:\n  {line.strip()}\n"
                    "Check that in_fmt matches the --outputcmd used in DetectNodes."
                )

            for offset, name in enumerate(columns):
                current[name].append(float(parts[n_leading + offset]))

            year, month, day, hour = (int(v) for v in parts[-n_date:])
            current["time"].append(datetime(year, month, day, hour))

    # Convert the accumulated lists to arrays once, at the end.
    for track in tracks:
        for name in columns:
            track[name] = np.asarray(track[name])
        track["time"] = np.asarray(track["time"])

    return tracks


def to_pm180(lon):
    """Convert longitude from the 0-360 convention TempestExtremes uses to
    -180..180, which is what matplotlib/cartopy expect for a CONUS-centred map.
    """
    return np.where(lon > 180.0, lon - 360.0, lon)


def split_at_dateline(lon, lat, threshold=180.0):
    """Yield (lon, lat) segments, broken wherever a track jumps the dateline.

    Without this, a storm crossing 180 degrees draws a spurious horizontal line
    straight across the whole map.
    """
    if len(lon) < 2:
        yield lon, lat
        return

    breaks = np.where(np.abs(np.diff(lon)) > threshold)[0] + 1
    start = 0
    for stop in breaks:
        yield lon[start:stop], lat[start:stop]
        start = stop
    yield lon[start:], lat[start:]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

# Map bounds, in -180..180. Wide enough to show the full life cycle of storms
# that pass over the US, not just their CONUS segment.
EXTENT = [-170.0, -30.0, 10.0, 75.0]

# The CONUS selection box applied in StitchNodes (--threshold), for reference.
CONUS_BOX = dict(lon=(234.0 - 360.0, 294.0 - 360.0), lat=(24.0, 52.0))


def _basemap(ax, title):
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="0.94")
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white")
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.4)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.3)
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.2, edgecolor="0.6")
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.7", alpha=0.6)
    gl.top_labels = False
    gl.right_labels = False
    ax.set_title(title, fontsize=11)


def _draw_conus_box(ax):
    lon0, lon1 = CONUS_BOX["lon"]
    lat0, lat1 = CONUS_BOX["lat"]
    ax.plot(
        [lon0, lon1, lon1, lon0, lon0],
        [lat0, lat0, lat1, lat1, lat0],
        transform=ccrs.PlateCarree(),
        color="tab:red", linewidth=1.2, linestyle="--", zorder=5,
    )


def plot_tracks(ax, tracks):
    """Panel 1: every storm track, coloured by its minimum central pressure."""
    _basemap(ax, f"ETC tracks (n = {len(tracks)})")
    _draw_conus_box(ax)

    norm = matplotlib.colors.Normalize(vmin=955.0, vmax=1010.0)
    cmap = plt.get_cmap("viridis_r")

    for track in tracks:
        lon = to_pm180(track["lon"])
        lat = track["lat"]

        # Colour the whole track by how deep the storm got (hPa).
        min_slp = np.min(track["slp"]) / 100.0
        color = cmap(norm(min_slp))

        for seg_lon, seg_lat in split_at_dateline(lon, lat):
            if len(seg_lon) > 1:
                ax.plot(seg_lon, seg_lat, transform=ccrs.PlateCarree(),
                        color=color, linewidth=1.0, alpha=0.85)

        # Mark genesis (first detected point) with a small dot.
        ax.plot(lon[0], lat[0], transform=ccrs.PlateCarree(),
                marker="o", markersize=2.5, color=color)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = plt.colorbar(sm, ax=ax, orientation="horizontal", pad=0.06, shrink=0.85)
    cb.set_label("Minimum central pressure (hPa)", fontsize=9)
    cb.ax.tick_params(labelsize=8)


def plot_density(ax, tracks, dlon=2.5, dlat=2.5):
    """Panel 2: track density -- how many 6-hourly storm points per grid box.

    This is the standard way ETC climatologies get compared between a model and
    a reanalysis, and it is what HistogramNodes computes for you. Doing it by
    hand here makes the definition explicit: it counts POINTS, so a slow-moving
    storm contributes more than a fast one.
    """
    _basemap(ax, "Track density (6-hourly points per grid box)")
    _draw_conus_box(ax)

    all_lon = to_pm180(np.concatenate([t["lon"] for t in tracks]))
    all_lat = np.concatenate([t["lat"] for t in tracks])

    lon_edges = np.arange(EXTENT[0], EXTENT[1] + dlon, dlon)
    lat_edges = np.arange(EXTENT[2], EXTENT[3] + dlat, dlat)

    counts, _, _ = np.histogram2d(all_lon, all_lat, bins=[lon_edges, lat_edges])
    counts = np.ma.masked_where(counts == 0, counts)

    mesh = ax.pcolormesh(
        lon_edges, lat_edges, counts.T,
        transform=ccrs.PlateCarree(), cmap="YlOrRd", shading="flat", zorder=3,
    )
    cb = plt.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.06, shrink=0.85)
    cb.set_label("Storm points per grid box", fontsize=9)
    cb.ax.tick_params(labelsize=8)


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("trajfile", help="StitchNodes trajectory file")
    parser.add_argument("--in_fmt", default="lon,lat,slp",
                        help="Column format, matching StitchNodes --in_fmt "
                             "(default: %(default)s)")
    parser.add_argument("--out", default="etc_tracks.png", help="Output PNG")
    args = parser.parse_args()

    tracks = read_trajectories(args.trajfile, in_fmt=args.in_fmt)

    if not tracks:
        raise SystemExit(
            f"No tracks found in {args.trajfile}. If DetectNodes ran but "
            "StitchNodes returned nothing, your --threshold or --mintime is "
            "probably too strict for the length of period you ran."
        )

    # A few numbers worth printing -- sanity checks before you trust any figure.
    lengths_h = [(t["time"][-1] - t["time"][0]).total_seconds() / 3600.0
                 for t in tracks]
    min_slps = [np.min(t["slp"]) / 100.0 for t in tracks]

    print(f"Storms            : {len(tracks)}")
    print(f"Mean duration     : {np.mean(lengths_h):.1f} h "
          f"(min {np.min(lengths_h):.0f}, max {np.max(lengths_h):.0f})")
    print(f"Mean min pressure : {np.mean(min_slps):.1f} hPa "
          f"(deepest {np.min(min_slps):.1f})")
    print(f"Total storm points: {sum(len(t['lon']) for t in tracks)}")

    fig, axes = plt.subplots(
        2, 1, figsize=(10, 11),
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    plot_tracks(axes[0], tracks)
    plot_density(axes[1], tracks)

    fig.suptitle(f"Extratropical cyclones: {args.trajfile}", fontsize=12)
    # Leave room on the left for the latitude labels, which cartopy draws
    # outside the axes and which bbox_inches="tight" will happily clip.
    fig.subplots_adjust(left=0.08, right=0.97, top=0.94, bottom=0.04, hspace=0.15)
    fig.savefig(args.out, dpi=140)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
