# twirl

twirl is a astrometric plate solving package for Python. It is suited for cases where the Right Ascension and Declination (RA, dec) coordinates of the image center and the field of view is known, computing a World Coordinate System (WCS) based on GAIA reference stars.

twirl is based on the algorithm of Lang et al. 2009 (astrometry.net) and go through these steps:

1. detection of stars in the image if not provided
2. catalog query using image known center
3. 4-points asterisms building and matching following Lang et al. 2009 method
4. image recombination and wcs fit using astropy.wcs

An offline version is under development (current version relies on a Gaia catalog query)

## Installation

```shell
pip install twirl
```

## Example
```python
import twirl

stars = twirl.find_peaks(data)[0:15]
wcs = twirl.compute_wcs(stars, center, fov)
```
A more complete example is provided in [docs/notebooks](https://github.com/lgrcia/twirl/tree/master/docs/notebooks)

## Development

### Project Requirements

- [Python](https://www.python.org/) 3.11.*
- [Poetry](https://python-poetry.org/) for Python package and environment management.

### Installing Dependencies

The twirl project manages Python package dependencies using [Poetry](https://python-poetry.org/). You'll need to follow the instructions for installation there.

Then you can start a shell session with the new environment with:

```console
$ poetry shell
```

**N.B.** For development with vscode you will need to run the following command:

```console
$ poetry config virtualenvs.in-project true
```

This will installed the poetry `.venv` in the root of the project and allow vscode to setup the environment correctly for development.

To start development, install all of the dependencies as:

```console
$ poetry install
```

**N.B.** _Ensure that any dependency changes are committed to source control, so everyone has a consistenct package dependecy list._

## Acknowledgements

This package has made use of Lang, D., Hogg, D.W., Mierle, K., Blanton, M. and Roweis, S. (2010). _Astrometry.net: Blind Astrometric Calibration of Arbitrary Astronomical Images_. The Astronomical Journal, 139(5), pp.1782–1800. [doi:10.1088/0004-6256/139/5/1782](https://iopscience.iop.org/article/10.1088/0004-6256/139/5/1782).