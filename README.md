# twirl

<p align="center">
    <img src="https://github.com/lgrcia/twirl/blob/main/docs/_static/twirl.png" height="270">
</p>

<p align="center">
  Python package for astrometric plate solving
  <br>
  <p align="center">
    <a href="https://github.com/lgrcia/twirl">
      <img src="https://img.shields.io/badge/github-lgrcia/twirl-blue.svg?style=flat" alt="github"/>
    </a>
    <a href="">
      <img src="https://img.shields.io/badge/license-MIT-lightgray.svg?style=flat" alt="license"/>
    </a>
    <a href="https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.4817G">
      <img src="https://img.shields.io/badge/paper-gray.svg?style=flat" alt="paper"/>
    </a>
    <a href="https://twirl.readthedocs.io">
      <img src="https://img.shields.io/badge/documentation-black.svg?style=flat" alt="documentation"/>
    </a>
  </p>
</p>

twirl is an astrometric plate solving package for Python. It is suited for cases where the Right Ascension and Declination (RA, dec) coordinates of the image center and the field of view is known, computing a World Coordinate System (WCS) based on GAIA reference stars.

twirl compute a WCS following these steps:

1. detection of stars in the image if not provided
2. catalog query using image known center
3. asterisms building and matching
4. image recombination and wcs fit using astropy.wcs

Astersisms are made of 3 or 4 points. 4 points asterisms are built following Lang et al. 2009 while 3 points asterims are based on an original algorithm.

## Installation

twirl can be installed using pip:

```shell
pip install twirl
```

or using poetry:

```shell
poetry add twirl
```

## Example Usage

*twirl* is designed to be complementary to the *astropy* package. It is used to compute a WCS by matching an image detected stars with catalog stars. To query the catalog stars, the center coordinates of the image must be known, as well as the size of the field of view. If not available, the latter can be computed using the known pixel scale of the detector.

Hence, the process starts by extracting the image RA-DEC center equatorial coordinate and compute the instrument field of view

```python
import numpy as np

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord

# Open some FITS image
hdu = fits.open("...")[0]

# get the center of the image
ra, dec = hdu.header["RA"], hdu.header["DEC"]
center = SkyCoord(ra, dec, unit=["deg", "deg"])

# and the size of its field of view
pixel = 0.66 * u.arcsec  # known pixel scale
shape = hdu.data.shape
fov = np.max(shape) * pixel.to(u.deg)
```

We can then query the gaia stars in the field using this information

```python
import twirl

sky_coords = twirl.gaia_radecs(center, 1.2 * fov)[0:12]
```

and match the queried stars to stars detected in the image

```python
# detect stars in the image
pixel_coords = twirl.find_peaks(hdu.data)[0:12]

# compute the World Coordinate System
wcs = twirl.compute_wcs(pixel_coords, sky_coords)
```
leading to a World Coordinate System object.

A more complete example is provided in [docs/ipynb/wcs.ipynb](https://twirl.readthedocs.io/en/latest/ipynb/wcs.html)


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

This package has made use of the algorithm from

Lang, D. et al. (2010). _Astrometry.net: Blind Astrometric Calibration of Arbitrary Astronomical Images_. The Astronomical Journal, 139(5), pp.1782–1800. [doi:10.1088/0004-6256/139/5/1782](https://iopscience.iop.org/article/10.1088/0004-6256/139/5/1782).

implemented in

Garcia, L. J. et al. (2022). prose: a Python framework for modular astronomical images processing. MNRAS, vol. 509, no. 4, pp. 4817–4828, 2022. [doi:10.1093/mnras/stab3113](https://academic.oup.com/mnras/article-abstract/509/4/4817/6414007).

See this [documentation page](https://twirl.readthedocs.io/en/latest/md/acknowledgement.html) for the BibTeX entries.
