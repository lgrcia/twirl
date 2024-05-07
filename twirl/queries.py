from datetime import datetime
from typing import Optional, Tuple, Union

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.units import Quantity


def gaia_radecs(
    center: Union[Tuple[float, float], SkyCoord],
    fov: Union[float, Quantity],
    limit: int = 10000,
    circular: bool = True,
    tmass: bool = False,
    dateobs: Optional[datetime] = None,
    magnitude: bool = False,
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
    tmass : bool, optional
        Whether to retrieve the 2MASS J magnitudes catelog. By default, it is set to False.
    dateobs : datetime.datetime, optional
        The date of the observation. If given, the proper motions of the sources will be taken into account. By default, it is set to None.

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

    fields = f"gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec {',gaia.phot_g_mean_mag' if magnitude else ''}"

    if circular and not tmass:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} {fields}
            FROM gaiadr2.gaia_source AS gaia
            WHERE 1=CONTAINS(
                POINT('ICRS', {ra}, {dec}), 
                CIRCLE('ICRS', gaia.ra, gaia.dec, {radius}))
            ORDER BY gaia.phot_g_mean_mag
            """
        )
    elif circular and tmass:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} {fields}
            FROM gaiadr2.gaia_source AS gaia
            INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match ON tmass_match.source_id = gaia.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = tmass_match.tmass_oid
            WHERE 1=CONTAINS(
                POINT('ICRS', {ra}, {dec}), 
                CIRCLE('ICRS', gaia.ra, gaia.dec, {radius}))
            ORDER BY tmass.j_m
            """
        )
    elif not circular and tmass:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} {fields}
            FROM gaiadr2.gaia_source AS gaia
            INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match ON tmass_match.source_id = gaia.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = tmass_match.tmass_oid
            WHERE gaia.ra BETWEEN {ra-ra_fov/2} AND {ra+ra_fov/2} AND
            gaia.dec BETWEEN {dec-dec_fov/2} AND {dec+dec_fov/2}
            ORDER BY tmass.j_m
            """
        )
    else:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec
            FROM gaiadr2.gaia_source AS gaia
            WHERE gaia.ra BETWEEN {ra-ra_fov/2} AND {ra+ra_fov/2} AND
            gaia.dec BETWEEN {dec-dec_fov/2} AND {dec+dec_fov/2}
            ORDER BY gaia.phot_g_mean_mag
            """
        )

    table = job.get_results()

    # add proper motion to ra and dec
    if dateobs is not None:
        # calculate fractional year
        dateobs = dateobs.year + (dateobs.timetuple().tm_yday - 1) / 365.25  # type: ignore

        years = dateobs - 2015.5  # type: ignore
        table["ra"] += years * table["pmra"] / 1000 / 3600
        table["dec"] += years * table["pmdec"] / 1000 / 3600

    radecs = np.array([table["ra"].value.data, table["dec"].value.data]).T

    if magnitude:
        return radecs, table["phot_g_mean_mag"].value.data
    else:
        return radecs
