from .utils import match
import numpy as np
from astropy.wcs.utils import fit_wcs_from_points
from astropy.coordinates import SkyCoord
import astropy.units as u


def gaia_radec(coord, shape, pixel, n=2000):
    """
    given the SkyCoord of target (astropy), the shape of image and the pixel size, return n gaia radec orderded by mag
    """
    radius = np.sqrt(2) * np.max(shape) * pixel / 120

    # querying gaia aroud target
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = n
    radius = u.Quantity(radius, u.arcminute)
    gaia_query = Gaia.cone_search_async(coord, radius, verbose=False)
    gaia_query_result = gaia_query.get_results()

    # ra dec result
    radec_gaia = np.array([
        gaia_query_result["ra"].data.data,
        gaia_query_result["dec"].data.data])

    # sorting in magnitudes
    idxs = np.argsort(gaia_query_result["phot_g_mean_flux"].data.data)[::-1]
    return radec_gaia[:, idxs].T


def match_wcs(xy, radec, return_aligned=False):
    """
    from list of pixel coords and radec coords performs the full wcs fit
    """
    matches, aligned = match(xy, radec, tolerance=20, return_aligned=True)
    radec = SkyCoord(*radec[matches[:, 0]].T, unit="deg")
    wcs = fit_wcs_from_points(xy[matches[:, 1]].T, radec)

    if return_aligned:
        return wcs, aligned
    else:
        return wcs


def find_peaks(data, threshold=2):
    from skimage.measure import label, regionprops
    threshold = threshold * np.nanstd(data.flatten()) + np.median(data.flatten())
    regions = regionprops(label(data > threshold), data)
    coordinates = np.array([region.weighted_centroid[::-1] for region in regions])
    fluxes = np.array([np.sum(region.intensity_image) for region in regions])

    return coordinates[np.argsort(fluxes)[::-1]]