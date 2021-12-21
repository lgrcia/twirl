# twirl

twirl is a python-only astrometric plate solving package. It is suited for cases where the RA-DEC coordinates of the image center and the field of view is known, computing a WCS based on gaia reference stars.

twirl is based on the algorithm of Lang et al. 2009 (astrometry.net) and go through these steps:
1. detection of stars in the image if not provided
2. catalog query using image known center
3. 4-points asterisms building and matching following Lang et al. 2009
4. image recombination and wcs fit (astropy.wcs)


An offline version is under development (current version relies on a Gaia catalog query)

## Example
```python
import twirl

stars = twirl.find_peaks(data)[0:15]
wcs = twirl.compute_wcs(stars, center, fov)
```
A more complete example is provided in [docs/notebooks](https://github.com/lgrcia/twirl/tree/master/docs/notebooks)

## Installation

```shell
pip install twirl
```