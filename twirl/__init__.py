from twirl.match import find_transform, get_transform_matrix, cross_match
from twirl.geometry import pad
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import fit_wcs_from_points
import numpy as np
import astropy.units as u

def compute_wcs(xy: np.ndarray, radecs: np.ndarray, tolerance: int=5):
    """Compute WCS based on unordered pixel vs. sky coordinates

    Parameters
    ----------
    xy : np.ndarray
        pixel coordinates
    radecs : np.ndarray
        RA-DEC coordinates (in deg)
    tolerance : int, optional
        minimum distance (in units of xy) between points to be cross-matched, by default 5

    Returns
    -------
    astropy.wcs.WCS
        image WCS 
    """
    M = find_transform(radecs, xy, tolerance=tolerance)
    radecs_xy = (M @ pad(radecs).T)[0:2].T
    i, j = cross_match(xy, radecs_xy).T
    M = get_transform_matrix(radecs[j], xy[i])
    radecs_xy = (M @ pad(radecs).T)[0:2].T
    i, j = cross_match(xy, radecs_xy).T
    return fit_wcs_from_points(xy[i].T, SkyCoord(radecs[j], unit="deg"))

def gaia_radecs(center, fov, limit=10000, circular=True):
    """
    Return RA-DEC (deg) in image based in image center and field-of-view

    
    query from https://gea.esac.esa.int/archive/documentation/GEDR3/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html

    TODO: adapt to non-square images

    Parameters
    ----------
    center : tuple or astropy.coordinates.SkyCoord
        image center sky coordinates (deg)
    fov : float or astropy.units.Unit
        field-of-view (deg)
    limit : int, optional
        limit queried sources, by default 10000
    circular : bool, optional
        whether query is circular, by default True

    Returns
    -------
    np.ndarray
        RA-DEC, shape (n, 2)
    """
    
    from astroquery.gaia import Gaia
    
    if isinstance(center, SkyCoord):
        ra = center.ra.to(u.deg).value
        dec = center.dec.to(u.deg).value
    
    if not isinstance(fov, u.Quantity):
        fov = fov * u.deg
    
    if fov.ndim == 1:
        ra_fov, dec_fov = fov.to(u.deg).value
    else:
        ra_fov = dec_fov = fov.to(u.deg).value

    radius = np.min([ra_fov, dec_fov])/2

    fields = 'ra, dec'

    if circular:
        job = Gaia.launch_job(f"select top {limit} {fields} from gaiadr2.gaia_source where "
                            "1=CONTAINS("
                            f"POINT('ICRS', {ra}, {dec}), "
                            f"CIRCLE('ICRS',ra, dec, {radius}))"
                            "order by phot_g_mean_mag")
    else:
        job = Gaia.launch_job(f"select top {limit} {fields} from gaiadr2.gaia_source where "
                    f"ra BETWEEN {ra-ra_fov/2} AND {ra+ra_fov/2} AND "
                    f"dec BETWEEN {dec-dec_fov/2} AND {dec+dec_fov/2} "
                    "order by phot_g_mean_mag") 

    table = job.get_results()
    return np.array([table["ra"].value.data, table["dec"].value.data]).T