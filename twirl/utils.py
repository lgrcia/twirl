import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astropy.wcs.utils import WCS, fit_wcs_from_points
from scipy.ndimage import gaussian_filter
from skimage.measure import label, regionprops

from twirl.geometry import pad
from twirl.match import cross_match, find_transform, get_transform_matrix
from twirl.queries import gaia_radecs


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
