"""
Measure how plate solving scales with the number of stars, and plot it.

For n stars in both lists (catalog and detections, all detected, no spurious
detections), this measures:

- the number of asterism hashes and of candidate hash pairs that
  `find_transform` must test
- the solve time of `compute_wcs`, with early stop (`min_match=0.8`) and
  without it (`min_match=None`, every pair is tested: the worst case)
- the cost of one tested pair (the body of the NumPy `find_transform` loop,
  timed on a sample of pairs)

The first two panels use log-log axes, so a power law n^k is a straight line
of slope k.

Usage
-----
    uv run python benchmarks/complexity.py
    uv run python benchmarks/complexity.py --output complexity.png --save data.json
    uv run python benchmarks/complexity.py --load data.json
    uv run python benchmarks/complexity.py --backend numpy
"""

import argparse
import json
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from scipy.spatial import cKDTree

from plate_solve import Case, make_field
from twirl import compute_wcs
from twirl.geometry import get_transform_matrix, pad
from twirl.match import count_cross_match
from twirl.quads import hashes as hash4
from twirl.triangles import hashes as hash3
from twirl.utils import _project_tangent_plane

ASTERISMS = {
    4: dict(label="quads", hash=hash4, color="#2a78d6"),
    3: dict(label="triangles", hash=hash3, color="#eb6834"),
}

# text and grid colors
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"


def pair_work(pixels, radecs, asterism, quads_tolerance, tolerance, sample=1000):
    """
    Hashes and candidate pairs as computed in `find_transform`, and the mean
    time (seconds) of one loop iteration, timed on up to `sample` pairs.
    """
    center = SkyCoord(*radecs.mean(0), unit="deg")
    radecs = _project_tangent_plane(center, SkyCoord(radecs, unit="deg")).T
    hash_function = ASTERISMS[asterism]["hash"]

    hashes_pixels, asterism_pixels = hash_function(pixels)
    hashes_radecs, asterism_radecs = hash_function(radecs)
    ball_query = cKDTree(hashes_pixels).query_ball_tree(
        cKDTree(hashes_radecs), r=quads_tolerance
    )
    n_pairs = sum(len(j) for j in ball_query)

    pairs = [(i, k) for i, j in enumerate(ball_query) for k in j[:10]][:sample]
    t0 = time.perf_counter()
    for i, j in pairs:
        M = get_transform_matrix(asterism_radecs[j], asterism_pixels[i])
        test = (M @ pad(radecs).T)[0:2].T
        count_cross_match(pixels, test, tolerance)
    per_pair = (time.perf_counter() - t0) / max(len(pairs), 1)

    return len(hashes_pixels), n_pairs, per_pair


def timed_solve(pixels, radecs, **kwargs):
    t0 = time.perf_counter()
    compute_wcs(pixels, radecs, **kwargs)
    return time.perf_counter() - t0


def measure(
    asterism, ns, max_solve, seeds, quads_tolerance, tolerance=5, backend="auto"
):
    rows = []
    for n in ns:
        case = Case("sweep", n, asterism=asterism, detected=1.0, n_spurious=0)
        solve = n <= max_solve
        values = {
            k: [] for k in ("hashes", "pairs", "per_pair_us", "early_s", "full_s")
        }

        for seed in seeds:
            pixels, radecs, *_ = make_field(case, seed)
            h, p, t = pair_work(pixels, radecs, asterism, quads_tolerance, tolerance)
            values["hashes"].append(h)
            values["pairs"].append(p)
            values["per_pair_us"].append(1e6 * t)
            if solve:
                kwargs = dict(
                    asterism=asterism,
                    tolerance=tolerance,
                    quads_tolerance=quads_tolerance,
                    backend=backend,
                )
                values["early_s"].append(
                    timed_solve(pixels, radecs, min_match=0.8, **kwargs)
                )
                values["full_s"].append(
                    timed_solve(pixels, radecs, min_match=None, **kwargs)
                )

        row = {"n": n, **{k: float(np.median(v)) for k, v in values.items() if v}}
        rows.append(row)

        line = (
            f"{ASTERISMS[asterism]['label']:<9} n={n:3d} hashes={row['hashes']:8.0f} "
            f"pairs={row['pairs']:10.0f} per pair={row['per_pair_us']:5.1f}us"
        )
        if solve:
            line += f" early={row['early_s']:7.3f}s full={row['full_s']:7.3f}s"
        print(line, flush=True)
    return rows


def series(rows, key):
    rows = [r for r in rows if key in r]
    return np.array([r["n"] for r in rows]), np.array([r[key] for r in rows])


def slope(n, y, last=4):
    """Power-law exponent fitted on the `last` largest n."""
    return np.polyfit(np.log(n[-last:]), np.log(y[-last:]), 1)[0]


def style(ax, title, ylabel, logy=True):
    ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", color=INK, fontsize=11, fontweight="bold")
    ax.set_xlabel("number of stars n", color=INK_2)
    ax.set_ylabel(ylabel, color=INK_2)
    ax.grid(True, which="major", color=GRID, linewidth=0.8)
    ax.tick_params(colors=INK_2, which="both")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)
    ticks = [5, 10, 20, 30, 50]
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.xaxis.set_minor_locator(plt.NullLocator())


def plot_line(ax, n, y, color, label, dashed=False, slope_label=True):
    ax.plot(
        n,
        y,
        color=color,
        linewidth=2,
        marker="o",
        markersize=4,
        linestyle="--" if dashed else "-",
        label=label,
    )
    if slope_label:
        ax.annotate(
            f"$n^{{{slope(n, y):.1f}}}$",
            (n[-1], y[-1]),
            xytext=(6, 0),
            textcoords="offset points",
            va="center",
            fontsize=10,
            color=INK,
        )


def plot(data, output):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), facecolor=SURFACE)

    ax = axes[0]
    for asterism, rows in data.items():
        a = ASTERISMS[asterism]
        plot_line(
            ax, *series(rows, "pairs"), a["color"], f"{a['label']}: candidate pairs"
        )
        plot_line(
            ax,
            *series(rows, "hashes"),
            a["color"],
            f"{a['label']}: hashes",
            dashed=True,
        )
    style(ax, "Work: hashes and pairs to test", "count")

    ax = axes[1]
    for asterism, rows in data.items():
        a = ASTERISMS[asterism]
        plot_line(
            ax,
            *series(rows, "full_s"),
            a["color"],
            f"{a['label']}: every pair (worst case)",
        )
        plot_line(
            ax,
            *series(rows, "early_s"),
            a["color"],
            f"{a['label']}: early stop",
            dashed=True,
            slope_label=False,
        )
    ax.axhline(10, color=INK_2, linewidth=1, linestyle=":")
    style(ax, "Solve time (compute_wcs)", "seconds")
    ax.annotate(
        "10 s",
        (ax.get_xlim()[0], 10),
        xytext=(4, 4),
        textcoords="offset points",
        fontsize=9,
        color=INK_2,
    )

    ax = axes[2]
    for asterism, rows in data.items():
        a = ASTERISMS[asterism]
        plot_line(
            ax, *series(rows, "per_pair_us"), a["color"], a["label"], slope_label=False
        )
    style(ax, "Cost of one tested pair", "microseconds", logy=False)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.1)

    # legends in empty corners: below the 10 s line in the solve time panel
    for ax, anchor in zip(axes, ((1, 0), (0, 0.82), (1, 0))):
        loc = "upper left" if anchor == (0, 0.82) else "lower right"
        ax.legend(
            frameon=False, fontsize=8.5, labelcolor=INK, loc=loc, bbox_to_anchor=anchor
        )
        # room for the slope labels
        lo, hi = ax.get_xlim()
        ax.set_xlim(lo, hi * 1.25)

    fig.suptitle(
        "twirl plate solving: how cost grows with the number of stars n "
        "(log-log: slope = exponent)",
        x=0.01,
        ha="left",
        color=INK,
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=150, facecolor=SURFACE)
    print(f"saved plot to {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", default="complexity.png")
    parser.add_argument("--save", help="write measurements to this JSON file")
    parser.add_argument("--load", help="plot measurements from a --save JSON file")
    parser.add_argument("--seeds", type=int, default=3, help="fields per point")
    parser.add_argument("--quads-tolerance", type=float, default=0.1)
    parser.add_argument("--backend", default="auto", choices=["auto", "numba", "numpy"])
    parser.add_argument(
        "--max-quads", type=int, default=30, help="largest n solved with quads"
    )
    parser.add_argument(
        "--max-triangles", type=int, default=40, help="largest n solved with triangles"
    )
    args = parser.parse_args()

    if args.load:
        with open(args.load) as f:
            data = {int(k): v for k, v in json.load(f).items()}
        plot(data, args.output)
        return

    seeds = range(args.seeds)
    sweeps = {
        4: ([6, 8, 10, 12, 15, 18, 20, 22, 25, 28, 30, 35, 40], args.max_quads),
        3: ([6, 8, 10, 12, 15, 18, 20, 25, 30, 35, 40, 50, 60], args.max_triangles),
    }

    data = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for asterism, (ns, max_solve) in sweeps.items():
            data[asterism] = measure(
                asterism,
                ns,
                max_solve,
                seeds,
                args.quads_tolerance,
                backend=args.backend,
            )

    if args.save:
        with open(args.save, "w") as f:
            json.dump({str(k): v for k, v in data.items()}, f, indent=2)
        print(f"saved data to {args.save}")

    plot(data, args.output)


if __name__ == "__main__":
    main()
