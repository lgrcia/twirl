> in dev, not functional

# twirl

twirl is a python-only astrometric plate solving package. It is suited for cases where the RA-DEC coordinate of the image center is known and a WCS needs to be computed easily.

twirl is based on the algorithm of Lang et al. 2009 (astrometry.net) and go through these steps:
1. detection of stars in the image if not provided
2. catalog query using image known center
3. 4-points asterisms building and matching following Lang et al. 2009
4. image recombination and wcs fit (astropy.wcs)

## Example
```python
import twirl

stars = twirl.find_peaks(data)[0:15]
radec = twirl.gaia_radec(target, shape, pixel)
wcs, aligned = twirl.match_wcs(stars, radec, return_aligned=True)
```
More examples are provided in [docs/notebooks](https://github.com/lgrcia/twirl/tree/master/docs/notebooks)

## Installation

For now twirl needs to be installed locally (soon pip installable)
```shell
git clone https://github.com/lgrcia/twirl.git
python3 -m pip install -e twirl
```