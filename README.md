# twirl

<p align="center">
    <img src=" docs/_static/twirl.png" height="270">
</p>

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

## Example Usage

twirl is designed to be complementary to the astropy package. It is used to compute a WCS from a set of stars detected in an image. 

As a prerequisite, star detetction and plate solving is suited for when the image center and field of view are known. 

In this case, the image center and field of view can be provided as a SkyCoord object and a Quantity object respectively.

Use any specified header that has been stored on the FITS primary HDU to obtain the center equatorial coordinate and field of view, for this example we will assum "RA" and "DEC" are the keywords for the center equatorial coordinate.

### Setup

```python
import numpy as np

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord

# Open some FITS image: 
hdul = fits.open("...")

# ra, dec in degrees
ra, dec = header["RA"], header["DEC"]

# Provide the center as a SkyCoord object:
center = SkyCoord(ra, dec, unit=["deg", "deg"])

center = [center.ra.value, center.dec.value]

# Utilise the image shape and pixel size in arcseconds " to obtain the field of view in degrees:
shape = data.shape

# Pixel size in arcseconds:
pixel = 0.66 * u.arcsec

# Field of view in degrees:
fov = np.max(shape)*pixel.to(u.deg).value
```

From here, we can pass the data, the center equatorial coordinate and the field-of-view to twirl to compute the World Coordinate System (WCS):

### Twirl Usage

```python
import twirl

# Find some starts in the image:
stars = twirl.find_peaks(data)[0:15]

# Compute the World Coordinate System:
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