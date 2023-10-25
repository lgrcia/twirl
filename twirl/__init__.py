from typing import Tuple, Union

import astropy
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astropy.wcs.utils import WCS, fit_wcs_from_points
from skimage.measure import label, regionprops

from twirl.geometry import pad
from twirl.match import cross_match, find_transform, get_transform_matrix


def _project_tangent_plane(center, skycoords):
    # Calculate the tangent plane projection
    tangent_plane_coords = skycoords.transform_to(center.skyoffset_frame())

    # Access the projected coordinates
    tangent_ra = tangent_plane_coords.lon
    tangent_dec = tangent_plane_coords.lat

    return np.array([tangent_ra.deg, tangent_dec.deg])


def compute_wcs(
    pixel_coords: np.ndarray,
    radecs: np.ndarray,
    tolerance: int = 5,
    quads_tolerance: float = 0.1,
    asterism=4,
    min_match=0.8,
) -> WCS:
    """
    Compute the WCS solution for an image given pixel coordinates and some unordered RA-DEC values.

    Parameters
    ----------
    pixel_coords : np.ndarray
        Pixel coordinates of the sources in the image, shape (n, 2)
    radecs : np.ndarray
        RA-DEC coordinates of the sources in the image, shape (m, 2)
    tolerance : int, optional
        Tolerance for the matching algorithm, by default 5
    asterism : int, optional
        Number of sources to use for matching, by default 4
    min_match : int, optional
        Minimum number of matches required, by default None

    Returns
    -------
    astropy.wcs.WCS
        WCS solution for the image if a match can be computed, None otherwise.
        A match is considered to be computed if at least one source and one target
        star are located less than `tolerance` pixels away from each other.
    """
    original_radecs = radecs.copy()
    center = SkyCoord(*radecs.mean(0), unit="deg")
    radecs = _project_tangent_plane(center, SkyCoord(radecs, unit="deg")).T

    M = find_transform(
        radecs,
        pixel_coords,
        tolerance=tolerance,
        asterism=asterism,
        min_match=min_match,
        quads_tolerance=quads_tolerance,
    )
    if M is None:
        return None
    else:
        radecs_xy = (M @ pad(radecs).T)[0:2].T
        i, j = cross_match(pixel_coords, radecs_xy).T
        M = get_transform_matrix(radecs[j], pixel_coords[i])
        radecs_xy = (M @ pad(radecs).T)[0:2].T
        i, j = cross_match(pixel_coords, radecs_xy).T
        return fit_wcs_from_points(
            pixel_coords[i].T, SkyCoord(original_radecs[j], unit="deg")
        )


def gaia_radecs(
    center: Union[Tuple[float, float], SkyCoord],
    fov: Union[float, Quantity],
    limit: int = 10000,
    circular: bool = True,
) -> np.ndarray:
    """
    Query the Gaia archive to retrieve the RA-DEC coordinates of stars within a given field-of-view (FOV) centered on a given sky position.

    Parameters
    ----------
    center : tuple or astropy.coordinates.SkyCoord
        The sky coordinates of the center of the FOV. If a tuple is given, it should contain the RA and DEC in degrees.
    fov : float or astropy.units.Quantity
        The field-of-view of the FOV in degrees. If a float is given, it is assumed to be in degrees.
    limit : int, optional
        The maximum number of sources to retrieve from the Gaia archive. By default, it is set to 10000.
    circular : bool, optional
        Whether to perform a circular or a rectangular query. By default, it is set to True.

    Returns
    -------
    np.ndarray
        An array of shape (n, 2) containing the RA-DEC coordinates of the retrieved sources in degrees.

    Raises
    ------
    ImportError
        If the astroquery package is not installed.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> from twirl import gaia_radecs
    >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
    >>> fov = 0.1
    >>> radecs = gaia_radecs(center, fov)
    """
    from astroquery.gaia import Gaia

    if isinstance(center, SkyCoord):
        ra = center.ra.deg
        dec = center.dec.deg
    else:
        ra, dec = center

    if not isinstance(fov, u.Quantity):
        fov = fov * u.deg

    if fov.ndim == 1:
        ra_fov, dec_fov = fov.to(u.deg).value
    else:
        ra_fov = dec_fov = fov.to(u.deg).value

    radius = np.min([ra_fov, dec_fov]) / 2

    fields = "ra, dec"

    if circular:
        job = Gaia.launch_job(
            f"select top {limit} {fields} from gaiadr2.gaia_source where "
            "1=CONTAINS("
            f"POINT('ICRS', {ra}, {dec}), "
            f"CIRCLE('ICRS',ra, dec, {radius}))"
            "order by phot_g_mean_mag"
        )
    else:
        job = Gaia.launch_job(
            f"select top {limit} {fields} from gaiadr2.gaia_source where "
            f"ra BETWEEN {ra-ra_fov/2} AND {ra+ra_fov/2} AND "
            f"dec BETWEEN {dec-dec_fov/2} AND {dec+dec_fov/2} "
            "order by phot_g_mean_mag"
        )

    table = job.get_results()
    return np.array([table["ra"].value.data, table["dec"].value.data]).T


def find_peaks(data: np.ndarray, threshold: float = 2.0) -> np.ndarray:
    """
    Find the coordinates of the peaks in a 2D array.

    Parameters
    ----------
    data : np.ndarray
        The 2D array to search for peaks.
    threshold : float, optional
        The threshold (in unit of image standard deviation) above which a pixel is considered
        part of a peak, i.e.
        The default is 2.0.

    Returns
    -------
    np.ndarray
        An array of shape (N, 2) containing the (x, y) coordinates of the N peaks
        found in the input array. The peaks are sorted by decreasing flux.
    """
    threshold = threshold * np.nanstd(data) + np.nanmedian(data)
    regions = regionprops(label(data > threshold), data)
    coordinates = np.array([region.weighted_centroid[::-1] for region in regions])
    fluxes = np.array([np.sum(region.intensity_image) for region in regions])

    return coordinates[np.argsort(fluxes)[::-1]]
