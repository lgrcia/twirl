"""
Benchmark the speed of plate solving on synthetic star fields.

Each star field is made from a known (true) TAN WCS:

1. catalog stars are drawn at random RA-DEC positions inside the field of view
2. a fraction of them is "detected": projected to pixels with the true WCS,
   with gaussian position noise
3. spurious detections (not in the catalog) are added and the detections are
   shuffled

`twirl.compute_wcs` then solves each field. A solve is a success when the
solved WCS places the detected catalog stars within `--max-error` pixels of
their true (noise-free) positions.

No network access is needed (no Gaia query).

Usage
-----
    uv run python benchmarks/plate_solve.py
    uv run python benchmarks/plate_solve.py --quick
    uv run python benchmarks/plate_solve.py --save baseline.json
    uv run python benchmarks/plate_solve.py --compare baseline.json
    uv run python benchmarks/plate_solve.py --backend numpy
"""

import argparse
import json
import platform
import time
import warnings
from dataclasses import asdict, dataclass

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import twirl
from twirl import compute_wcs

IMAGE_SIZE = 2048  # pixels
PIXEL_SCALE = 0.66 / 3600  # degrees per pixel


@dataclass(frozen=True)
class Case:
    name: str
    n_stars: int  # catalog stars in the field
    asterism: int = 4
    detected: float = 0.8  # fraction of catalog stars that are detected
    n_spurious: int = 2  # detections that are not in the catalog
    noise: float = 0.5  # detection position noise (pixels, 1 sigma)


SUITE = [
    Case("perfect-12-quads", 12, asterism=4, detected=1.0, n_spurious=0, noise=0.0),
    Case("perfect-12-triangles", 12, asterism=3, detected=1.0, n_spurious=0, noise=0.0),
    Case("10-quads", 10, asterism=4),
    Case("10-triangles", 10, asterism=3),
    Case("15-quads", 15, asterism=4),
    Case("15-triangles", 15, asterism=3),
    Case("20-quads", 20, asterism=4),
    Case("20-triangles", 20, asterism=3),
    Case("30-quads", 30, asterism=4, n_spurious=4),
    Case("30-triangles", 30, asterism=3, n_spurious=4),
    Case("crowded-20-quads", 20, asterism=4, detected=0.6, n_spurious=8, noise=1.0),
    Case("crowded-20-triangles", 20, asterism=3, detected=0.6, n_spurious=8, noise=1.0),
]

QUICK_SUITE = [c for c in SUITE if c.n_stars <= 15]


def true_wcs(rng):
    """A TAN WCS with random sky center and random rotation."""
    ra = rng.uniform(0, 360)
    dec = np.rad2deg(np.arcsin(rng.uniform(-0.98, 0.98)))
    angle = rng.uniform(0, 2 * np.pi)
    c, s = np.cos(angle), np.sin(angle)

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.crpix = [IMAGE_SIZE / 2, IMAGE_SIZE / 2]
    # flip RA axis (east left) as on a real sky image
    wcs.wcs.cd = PIXEL_SCALE * np.array([[-c, s], [s, c]])
    return wcs


def make_field(case, seed):
    """
    Returns
    -------
    pixels : np.ndarray
        detected positions (pixels), shape (n_detections, 2)
    radecs : np.ndarray
        catalog positions (degrees), shape (n_stars, 2)
    wcs : astropy.wcs.WCS
        the true WCS
    detected_radecs, detected_pixels : np.ndarray
        catalog stars that are detected and their noise-free pixel positions
    """
    rng = np.random.default_rng(seed)
    wcs = true_wcs(rng)

    # catalog stars, drawn uniformly in the image area
    xy = rng.uniform(0, IMAGE_SIZE, size=(case.n_stars, 2))
    radecs = np.array(wcs.pixel_to_world_values(*xy.T)).T

    n_detected = max(4, int(round(case.detected * case.n_stars)))
    idx = rng.choice(case.n_stars, n_detected, replace=False)
    detected_pixels = xy[idx]
    noisy = detected_pixels + rng.normal(0, case.noise, size=detected_pixels.shape)
    spurious = rng.uniform(0, IMAGE_SIZE, size=(case.n_spurious, 2))
    pixels = np.vstack([noisy, spurious])
    rng.shuffle(pixels)

    return pixels, radecs, wcs, radecs[idx], detected_pixels


def solve_error(wcs, radecs, pixels):
    """Median distance (pixels) between where `wcs` puts `radecs` and `pixels`."""
    if wcs is None:
        return np.inf
    xy = np.array(wcs.world_to_pixel(SkyCoord(radecs, unit="deg"))).T
    return float(np.median(np.linalg.norm(xy - pixels, axis=1)))


def run_case(case, seeds, repeat, max_error, solver_kwargs):
    times = []  # best time of each field (seconds)
    errors = []

    for seed in seeds:
        pixels, radecs, _, det_radecs, det_pixels = make_field(case, seed)

        best = np.inf
        solved = None
        for _ in range(repeat):
            t0 = time.perf_counter()
            try:
                solved = compute_wcs(
                    pixels, radecs, asterism=case.asterism, **solver_kwargs
                )
            except Exception:
                solved = None
            best = min(best, time.perf_counter() - t0)

        times.append(best)
        errors.append(solve_error(solved, det_radecs, det_pixels))

    times = np.array(times)
    errors = np.array(errors)
    return {
        **asdict(case),
        "n_fields": len(seeds),
        "median_ms": 1e3 * float(np.median(times)),
        "min_ms": 1e3 * float(np.min(times)),
        "max_ms": 1e3 * float(np.max(times)),
        "total_ms": 1e3 * float(np.sum(times)),
        "success_rate": float(np.mean(errors < max_error)),
        "median_error_px": float(np.median(errors)),
    }


def print_table(results, baseline=None):
    header = (
        f"{'case':<22} {'stars':>5} {'det':>4} {'median ms':>10} "
        f"{'min ms':>8} {'max ms':>8} {'solved':>7}"
    )
    if baseline:
        header += f" {'base ms':>8} {'speedup':>8}"
    print(header)
    print("-" * len(header))

    for r in results:
        n_det = max(4, int(round(r["detected"] * r["n_stars"]))) + r["n_spurious"]
        line = (
            f"{r['name']:<22} {r['n_stars']:>5} {n_det:>4} {r['median_ms']:>10.2f} "
            f"{r['min_ms']:>8.2f} {r['max_ms']:>8.2f} {r['success_rate']:>7.0%}"
        )
        if baseline:
            base = baseline.get(r["name"])
            if base is None:
                line += f" {'-':>8} {'-':>8}"
            else:
                speedup = base["median_ms"] / r["median_ms"]
                line += f" {base['median_ms']:>8.2f} {speedup:>7.2f}x"
        print(line)

    total = sum(r["total_ms"] for r in results)
    print("-" * len(header))
    print(f"total solve time (best of repeats, all fields): {total / 1e3:.2f} s")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--fields", type=int, default=20, help="fields per case")
    parser.add_argument("--repeat", type=int, default=3, help="timings per field")
    parser.add_argument("--quick", action="store_true", help="small cases only")
    parser.add_argument("--case", action="append", help="run only these cases")
    parser.add_argument("--max-error", type=float, default=2.0, help="pixels")
    parser.add_argument("--tolerance", type=float, default=5)
    parser.add_argument("--quads-tolerance", type=float, default=0.1)
    parser.add_argument("--min-match", type=float, default=0.8)
    parser.add_argument("--save", help="write results to this JSON file")
    parser.add_argument("--compare", help="compare with a JSON file from --save")
    parser.add_argument("--backend", default="auto", choices=["auto", "numba", "numpy"])
    args = parser.parse_args()

    cases = QUICK_SUITE if args.quick else SUITE
    if args.case:
        cases = [c for c in SUITE if c.name in args.case]
        unknown = set(args.case) - {c.name for c in cases}
        if unknown:
            parser.error(f"unknown case(s): {', '.join(sorted(unknown))}")

    solver_kwargs = dict(
        tolerance=args.tolerance,
        quads_tolerance=args.quads_tolerance,
        min_match=args.min_match,
        backend=args.backend,
    )
    seeds = range(args.fields)

    print(
        f"twirl {twirl.__file__}\n"
        f"python {platform.python_version()}, numpy {np.__version__}, "
        f"{platform.machine()} {platform.processor()}\n"
        f"backend {args.backend}, "
        f"{args.fields} fields per case, best of {args.repeat} runs per field\n"
    )

    # warm up imports and caches so the first case is not penalized
    pixels, radecs, *_ = make_field(SUITE[0], 0)
    compute_wcs(pixels, radecs, backend=args.backend)

    results = []
    with warnings.catch_warnings():
        # fit_wcs_from_points and SkyCoord can warn on some fields
        warnings.simplefilter("ignore")
        for case in cases:
            results.append(
                run_case(case, seeds, args.repeat, args.max_error, solver_kwargs)
            )

    baseline = None
    if args.compare:
        with open(args.compare) as f:
            baseline = {r["name"]: r for r in json.load(f)["results"]}

    print_table(results, baseline)

    if args.save:
        with open(args.save, "w") as f:
            json.dump({"args": vars(args), "results": results}, f, indent=2)
        print(f"saved to {args.save}")


if __name__ == "__main__":
    main()
