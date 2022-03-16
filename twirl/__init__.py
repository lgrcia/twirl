from . import utils
import numpy as np
from astropy.wcs.utils import fit_wcs_from_points
from astropy.coordinates import SkyCoord
import astropy.units as u
from .config import ConfigManager

CONFIG = ConfigManager()

def gaia_radecs(center, fov, limit=3000):
    from astroquery.gaia import Gaia

    job = Gaia.launch_job(f"select top {limit} ra, dec from gaiadr2.gaia_source where "
                          "1=CONTAINS("
                          f"POINT('ICRS', {center[0]}, {center[1]}), "
                          f"CIRCLE('ICRS',ra, dec, {fov/2}))"
                          #f"ra between {center[0]-fov/2} and {center[0]+fov/2} and "
                          #f"dec between {center[1]-fov/2} and {center[1]+fov/2}"
                          "order by phot_g_mean_mag")

    table = job.get_results()
    return np.array([table["ra"].data.data, table["dec"].data.data]).T

def compute_wcs(stars, center, fov, offline=False, n=15):
    if offline:
        CONFIG.load_catalog()
        gaias = CONFIG.current_catalog
        isin = np.all(np.abs(gaias - center) < fov, 1)
        gaias = gaias[isin]
    else:
        gaias = gaia_radecs(center, fov)
    
    return _compute_wcs(stars, gaias, n=n)

def _compute_wcs(stars, gaias, n=15, tolerance=10):
    X = utils.find_transform(gaias[0:n], stars, n=n, tolerance=tolerance)
    gaia_pixels = utils.affine_transform(X)(gaias)
    s1, s2 = utils.cross_match(gaia_pixels, stars, return_ixds=True, tolerance=15).T
    ras, decs = gaias[s1].T
    gaia_coords = SkyCoord(ra=ras*u.deg, dec=decs*u.deg)
    return fit_wcs_from_points(stars[s2].T, gaia_coords)

def find_peaks(data, threshold=2):
    from skimage.measure import label, regionprops
    threshold = threshold * np.nanstd(data.flatten()) + np.median(data.flatten())
    regions = regionprops(label(data > threshold), data)
    coordinates = np.array([region.weighted_centroid[::-1] for region in regions])
    fluxes = np.array([np.sum(region.intensity_image) for region in regions])

    return coordinates[np.argsort(fluxes)[::-1]]


from pkg_resources import get_distribution
__version__ = get_distribution('twirl').version
