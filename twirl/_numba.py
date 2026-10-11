"""
Numba-compiled loop of `twirl.match.find_transform`.

This module is only imported when numba is installed and the "numba" backend
is used, so that `import twirl` stays fast and numba stays optional.
"""

import numba
import numpy as np


@numba.njit(cache=True)
def match_counts(
    asterism_radecs, asterism_pixels, i, j, radecs, pixels, tolerance, stop
):
    """
    For each candidate pair (i[p], j[p]), fit the affine transform from
    `asterism_radecs[j[p]]` to `asterism_pixels[i[p]]` (least squares), apply it
    to `radecs` and count the transformed points closer than `tolerance` to any
    of `pixels`.

    Same result as the loop body of `find_transform`. Stops after the first
    pair whose count is `>= stop` (if `stop >= 0`) and returns the counts up to
    and including that pair.
    """
    n_pairs = len(i)
    counts = np.zeros(n_pairs, np.int64)
    tolerance2 = tolerance * tolerance

    for p in range(n_pairs):
        src = asterism_radecs[j[p]]
        dst = asterism_pixels[i[p]]
        n = src.shape[0]

        # centered coordinates: the linear part is a 2x2 least squares problem
        # and the translation follows from the centroids
        mx = my = mu = mv = 0.0
        for q in range(n):
            mx += src[q, 0]
            my += src[q, 1]
            mu += dst[q, 0]
            mv += dst[q, 1]
        mx /= n
        my /= n
        mu /= n
        mv /= n

        sxx = sxy = syy = sxu = syu = sxv = syv = 0.0
        for q in range(n):
            x = src[q, 0] - mx
            y = src[q, 1] - my
            u = dst[q, 0] - mu
            v = dst[q, 1] - mv
            sxx += x * x
            sxy += x * y
            syy += y * y
            sxu += x * u
            syu += y * u
            sxv += x * v
            syv += y * v

        det = sxx * syy - sxy * sxy
        if det == 0.0:
            continue
        a = (syy * sxu - sxy * syu) / det
        b = (sxx * syu - sxy * sxu) / det
        c = (syy * sxv - sxy * syv) / det
        d = (sxx * syv - sxy * sxv) / det
        tx = mu - a * mx - b * my
        ty = mv - c * mx - d * my

        count = 0
        for r in range(radecs.shape[0]):
            x = a * radecs[r, 0] + b * radecs[r, 1] + tx
            y = c * radecs[r, 0] + d * radecs[r, 1] + ty
            for s in range(pixels.shape[0]):
                dx = pixels[s, 0] - x
                dy = pixels[s, 1] - y
                if dx * dx + dy * dy < tolerance2:
                    count += 1
                    break
        counts[p] = count

        if stop >= 0 and count >= stop:
            return counts[: p + 1]

    return counts
