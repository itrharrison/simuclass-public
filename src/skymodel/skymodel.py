"""
Script to convert a T-RECS catalogue into a sky model FITS file.
Forms part of the larger simuCLASS pipeline for creating simulated superCLASS
observations.

Usage:
python skymodel.py example.ini

Prequisities:
skymodel_tools.py
galsim
astropy

Contact:
Ian Harrison
ian.harrison-2@manchester.ac.uk

"""
import numpy as np
from numpy.core.defchararray import add as stradd
from numpy.core.defchararray import multiply as strmultiply
import pdb
import sys
import ConfigParser
import time
import cPickle as pickle

import galsim

from astropy.io import fits
from astropy.table import Table
from astropy import units as uns
from astropy import wcs as ast_wcs
from astropy.coordinates import SkyCoord

from scipy import random
from matplotlib import pyplot as plt

from skymodel_tools import setup_wcs

big_fft_params = galsim.GSParams(maximum_fft_size=84188)
arcsectorad = (1.0 * uns.arcsec).to(uns.rad).value
degtoarcsec = (1.0 * uns.deg).to(uns.arcsec).value


def runSkyModel(config):
    """Simulate a sky model from a T-RECS catalogue.

  Parameters
  ----------
  config : configparser
    ConfigParser configuration containing necessary sections.

  """
    output_path = config.get("pipeline", "output_dirname")

    # Set some image properties
    pixel_scale = config.getfloat("skymodel", "pixel_scale") * galsim.arcsec
    fov = config.getfloat("skymodel", "field_of_view") * galsim.arcmin
    image_size = int((fov / galsim.arcmin) / (pixel_scale / galsim.arcmin))

    ra_field = config.get("skymodel", "field_ra")
    ra_field_gs = galsim.HMS_Angle(ra_field)
    dec_field = config.get("skymodel", "field_dec")
    dec_field_gs = galsim.DMS_Angle(dec_field)

    ra_survey = config.get("skymodel", "survey_ra")
    ra_survey_gs = galsim.HMS_Angle(ra_survey)
    dec_survey = config.get("skymodel", "survey_dec")
    dec_survey_gs = galsim.DMS_Angle(dec_survey)

    w_twod = setup_wcs(config, ndim=2)
    w_fourd = setup_wcs(config, ndim=4)
    header_twod = w_twod.to_header()
    header_fourd = w_fourd.to_header()
    header_fourd["BUNIT"] = "JY/PIXEL"

    # Create the galsim image
    full_image = galsim.ImageF(image_size, image_size, scale=pixel_scale)
    im_center = full_image.bounds.trueCenter()
    sky_center = galsim.CelestialCoord(ra=ra_field_gs, dec=dec_field_gs)

    # Create a WCS for the galsim image
    full_image.wcs, origin = galsim.wcs.readFromFitsHeader(header_twod)

    tstart = time.time()

    if config.getboolean("skymodel", "dosfg"):
        # Load the catalogue
        cat_file_name = config.get("skymodel", "catalogue_filepath")
        print("Loading catalogue from {0} ...".format(cat_file_name))
        cat = Table()
        if config.get("skymodel", "type") == "trecs":
            cat_read = Table.read(cat_file_name, format="ascii")

            source_prefix = "TRECS-"
            source_name_ra = np.asarray(cat_read["lon"], dtype=str)
            source_name_dec = np.asarray(cat_read["lat"], dtype=str)

            source_prefix_arr = strmultiply(
                source_prefix, np.ones_like(source_name_ra, dtype=int)
            )
            source_l_arr = strmultiply("l", np.ones_like(source_name_ra, dtype=int))
            source_b_arr = strmultiply("b", np.ones_like(source_name_ra, dtype=int))

            source_name_pos = stradd(
                source_prefix_arr,
                stradd(
                    source_l_arr,
                    stradd(source_name_ra, (stradd(source_b_arr, source_name_dec))),
                ),
            )
            cat["Source_id"] = source_name_pos

            cat["ra_offset"] = cat_read["lon"]  # deg
            cat["ra_offset"].unit = "deg"

            cat["dec_offset"] = cat_read["lat"]  # deg
            cat["dec_offset"].unit = "deg"

            # convert the offsets from the *survey* centre to RA and DEC
            cat["DEC"] = dec_survey_gs / galsim.degrees + cat["dec_offset"]
            dec_abs_radians = cat["DEC"] * galsim.degrees / galsim.radians
            cat["RA"] = ra_survey_gs / galsim.degrees + cat["ra_offset"] / np.cos(
                np.asarray(dec_abs_radians, dtype=float)
            )

            # calculate the offsets from the centre of the *pointing*
            cat["dec_offset"] = cat["DEC"] - dec_field_gs / galsim.degrees
            dec_abs_radians = cat["DEC"] * galsim.degrees / galsim.radians
            cat["ra_offset"] = (cat["RA"] - ra_field_gs / galsim.degrees) * np.cos(
                np.asarray(dec_abs_radians, dtype=float)
            )

            cat["bulge_disk_amplitude_ratio"] = cat_read["bulge/disk"]

            cat["Total_flux"] = cat_read["flux"] * 1.0e-3  # Jy
            cat["Total_flux"].unit = "Jy"

            cat["Maj"] = cat_read["size"]  # arcsec
            cat["Maj"].unit = "arcsec"

            scale_radius_to_hlr = 1.6783469900166605
            cat["Maj_halflight"] = cat_read["size"] * scale_radius_to_hlr
            cat["Maj_halflight"].unit = "arcsec"

            cat["Peak_flux"] = cat["Total_flux"] / (2.0 * cat["Maj"] * arcsectorad)
            cat["Peak_flux"].unit = "Jy"

            cat["e1"] = cat_read["e1"]
            cat["e2"] = cat_read["e2"]

            cat["g1_shear"] = cat_read["gamma1"]
            cat["g2_shear"] = cat_read["gamma2"]

            cat["mod_e"] = np.sqrt(cat["e1"] ** 2.0 + cat["e2"] ** 2.0)

            cat["q"] = q_obs = (1.0 - cat["mod_e"] ** 2.0) / (1.0 + cat["mod_e"] ** 0.2)

            cat["Min"] = cat["Maj"] * cat["q"]
            cat["Min"].unit = "arcsec"

            cat["Min_halflight"] = cat["Maj_halflight"] * cat["q"]
            cat["Min_halflight"].unit = "arcsec"

            cat["PA"] = 0.5 * np.arctan2(cat["e2"], cat["e1"])
            cat["PA"].unit = "rad"

        elif config.get("skymodel", "type") == "pybdsm":
            cat_read = Table.read(cat_file_name, format="fits")

            cat["Source_id"] = cat_read["Source_id"]

            cat["RA"] = cat_read["RA"]  # deg
            cat["DEC"] = cat_read["DEC"]  # deg

            cat["dec_offset"] = cat["DEC"] - dec_field_gs / galsim.degrees
            dec_abs_radians = cat["DEC"] * galsim.degrees / galsim.radians
            cat["ra_offset"] = (cat["RA"] - ra_field_gs / galsim.degrees) * np.cos(
                np.asarray(dec_abs_radians, dtype=float)
            )

            cat["Total_flux"] = cat_read["Total_flux"]  # Jy
            cat["Total_flux"].unit = "Jy"

            cat["Maj"] = cat_read["Maj"] * degtoarcsec  # deg
            cat["Maj"].unit = "arcsec"

            cat["Maj_halflight"] = 0.5 * cat["Maj"]
            cat["Maj_halflight"].unit = "arcsec"

            cat["Min"] = cat_read["Min"] * degtoarcsec  # deg
            cat["Min"].unit = "arcsec"

            cat["Min_halflight"] = 0.5 * cat["Min"]
            cat["Min_halflight"].unit = "arcsec"

            cat["Peak_flux"] = cat_read["Peak_flux"]  # Jy
            cat["Peak_flux"].unit = "Jy"

            cat["q"] = cat_read["Min"] / cat_read["Maj"]
            cat["PA"] = (cat_read["PA"] * uns.deg).to(uns.rad).value
            cat["PA"].unit = "rad"
            cat["mod_e"] = (1.0 - cat["q"] ** 2.0) / (1.0 + cat["q"] ** 2.0)

            cat["e1"] = cat["mod_e"] * np.cos(2.0 * cat["PA"])
            cat["e2"] = cat["mod_e"] * np.sin(2.0 * cat["PA"])

            cat["g1_shear"] = 0.0e0
            cat["g2_shear"] = 0.0e0

        # fov cut
        ra_offset_max = 0.9 * (fov / 2) / galsim.degrees
        dec_offset_max = 0.9 * (fov / 2) / galsim.degrees

        fov_cut = (abs(cat["ra_offset"]) < ra_offset_max) * (
            abs(cat["dec_offset"]) < dec_offset_max
        )
        cat = cat[fov_cut]

        # flux cuts
        if config.getboolean("skymodel", "highfluxcut"):
            highflux_cut = cat["Total_flux"] < config.getfloat(
                "skymodel", "highfluxcut_value"
            )
            cat = cat[highflux_cut]

        if config.getboolean("skymodel", "lowfluxcut"):
            lowflux_cut = cat["Total_flux"] > config.getfloat(
                "skymodel", "lowfluxcut_value"
            )
            cat = cat[lowflux_cut]

        if config.getboolean("skymodel", "highsizecut"):
            highsize_cut = cat["Maj"] < config.getfloat("skymodel", "highsizecut_value")
            cat = cat[highsize_cut]

        if config.getboolean("skymodel", "lowsizecut"):
            lowsize_cut = cat["Maj"] > config.getfloat("skymodel", "lowsizecut_value")
            cat = cat[lowsize_cut]

        if config.get("skymodel", "sizescale") == "constant":
            cat["Maj"] = np.ones_like(cat["Maj"]) * config.getfloat(
                "skymodel", "sizescale_constant_value"
            )

            scale_radius_to_hlr = 1.6783469900166605
            cat["Maj_halflight"] = cat["Maj"] * scale_radius_to_hlr
            cat["Maj_halflight"].unit = "arcsec"

            cat["Min"] = cat["Maj"] * cat["q"]
            cat["Min"].unit = "arcsec"

            cat["Min_halflight"] = cat["Maj_halflight"] * cat["q"]
            cat["Min_halflight"].unit = "arcsec"

        # number of sources, on grid if requested
        if config.getboolean("skymodel", "grid"):
            nobj = int(np.sqrt(config.getint("skymodel", "ngals"))) ** 2.0
            cat["ra_offset"] = np.linspace(-ra_offset_max, ra_offset_max, nobj)
            cat["dec_offset"] = np.linspace(-ra_offset_max, ra_offset_max, nobj)
        else:
            nobj = len(cat)
            if config.getint("skymodel", "ngals") > -1:
                nobj = config.getint("skymodel", "ngals")
                cat = cat[:nobj]

        # flux range
        if config.get("skymodel", "fluxscale") == "constant":
            cat["Total_flux"] = np.ones_like(cat["Total_flux"]) * config.getfloat(
                "skymodel", "fluxscale_constant_value"
            )
            cat["Peak_flux"] = cat["Total_flux"] / (2.0 * cat["Maj"] * arcsectorad)

        # scale flux
        cat["Total_flux"] = cat["Total_flux"] * config.getfloat(
            "skymodel", "flux_factor"
        )
        cat["Peak_flux"] = cat["Peak_flux"] * config.getfloat("skymodel", "flux_factor")

        # scale size
        cat["Maj"] = cat["Maj"] * config.getfloat("skymodel", "sizefactor")
        cat["Maj_halflight"] = cat["Maj_halflight"] * config.getfloat(
            "skymodel", "sizefactor"
        )
        cat["Min"] = cat["Min"] * config.getfloat("skymodel", "sizefactor")
        cat["Min_halflight"] = cat["Min_halflight"] * config.getfloat(
            "skymodel", "sizefactor"
        )

        # write out catalogue

        output_filename = "truthcat" + config.get("pipeline", "output_suffix") + ".fits"

        print("Writing truthcat to " + os.path.join(output_path, output_filename))
        cat.write(
            os.path.join(output_path, output_filename), format="fits", overwrite=True,
        )

        ix_arr = np.ones(nobj)
        iy_arr = np.ones(nobj)
        print("...done.")

        # Draw the galaxies onto the galsim image
        for i, cat_gal in enumerate(cat):

            sys.stdout.write(
                "\rAdding source {0} of {1} to skymodel...".format(i + 1, nobj)
            )

            # choose the profile
            if config.get("skymodel", "galaxy_profile") == "exponential":
                gal = galsim.Exponential(
                    scale_radius=cat_gal["Maj"] / 2.0,
                    flux=cat_gal["Total_flux"],
                    gsparams=big_fft_params,
                )

            elif config.get("skymodel", "galaxy_profile") == "gaussian":
                gal = galsim.Gaussian(
                    fwhm=cat_gal["Maj"],
                    flux=cat_gal["Total_flux"],
                    gsparams=big_fft_params,
                )

            elif config.get("skymodel", "galaxy_profile") == "matched-exponential":
                gauss_gal = galsim.Gaussian(
                    fwhm=cat_gal["Maj"], flux=cat_gal["Total_flux"]
                )
                gal = galsim.Exponential(
                    half_light_radius=gauss_gal.getHalfLightRadius(),
                    flux=cat_gal["Total_flux"],
                    gsparams=big_fft_params,
                )
                print("scale radius: {}".format(gal.getScaleRadius()))
                del gauss_gal

            elif config.get("skymodel", "galaxy_profile") == "sersic":
                gauss_gal = galsim.Gaussian(
                    fwhm=cat_gal["Maj"], flux=cat_gal["Total_flux"]
                )
                gal = galsim.Sersic(
                    n=config.getfloat("skymodel", "sersic_index"),
                    half_light_radius=gauss_gal.getHalfLightRadius(),
                    flux=cat_gal["Total_flux"],
                    gsparams=big_fft_params,
                )
                del gauss_gal

            elif config.get("skymodel", "galaxy_profile") == "composite":
                gauss_gal = galsim.Gaussian(
                    fwhm=cat_gal["Maj"], flux=cat_gal["Total_flux"]
                )
                disk = galsim.Gaussian(
                    fwhm=cat_gal["Maj"],
                    flux=cat_gal["Total_flux"],
                    gsparams=big_fft_params,
                )
                bulge = galsim.Exponential(
                    half_light_radius=disk.getHalfLightRadius(),
                    flux=cat_gal["Total_flux"],
                    gsparams=big_fft_params,
                )
                gal = disk + cat_gal["bulge_disk_amplitude_ratio"] * bulge
                gal = gal.withFlux(cat_gal["Total_flux"])
                del gauss_gal
            elif config.get("skymodel", "galaxy_profile") == "points":
                gal = galsim.Gaussian(
                    fwhm=1.0e-2, flux=cat_gal["Total_flux"], gsparams=big_fft_params
                )

            if config.getboolean("skymodel", "truncate_profile"):
                gal = galsim.Sersic(
                    1,
                    scale_radius=cat_gal["Maj"] / 2,
                    flux=cat_gal["Total_flux"],
                    trunc=cat_gal["Maj"],
                    gsparams=big_fft_params,
                )

            # calculate the total ellipticity
            ellipticity = galsim.Shear(e1=cat_gal["e1"], e2=cat_gal["e2"])
            shear = galsim.Shear(g1=cat_gal["g1_shear"], g2=cat_gal["g2_shear"])
            if config.getboolean("skymodel", "doshear"):
                total_shear = ellipticity + shear
            else:
                total_shear = ellipticity

            maj_gal = cat_gal["Maj"]
            q_gal = cat_gal["Min"] / cat_gal["Maj"]
            A_gal = np.pi * maj_gal ** 2.0
            maj_corr_gal = np.sqrt(A_gal / (np.pi * q_gal))

            gal = gal.shear(total_shear)
            gal = gal.dilate(maj_gal / maj_corr_gal)

            if config.getboolean("skymodel", "dosimple_psf"):
                psf_maj = config.getfloat("skymodel", "simple_psf_maj") * galsim.arcsec
                psf_min = config.getfloat("skymodel", "simple_psf_min") * galsim.arcsec
                psf_pa = (
                    config.getfloat("skymodel", "simple_psf_pa") - 90.0
                ) * galsim.degrees
                q = (psf_min / galsim.arcsec) / (psf_maj / galsim.arcsec)
                psf = galsim.Gaussian(fwhm=psf_maj / galsim.arcsec)
                psf_shear = galsim.Shear(q=q, beta=psf_pa)
                psf = psf.shear(psf_shear)

                gal = galsim.Convolve(gal, psf)

            x, y = w_twod.wcs_world2pix(cat_gal["RA"], cat_gal["DEC"], 0,)
            x = float(x)
            y = float(y)

            # Account for the fractional part of the position:
            ix = int(np.floor(x + 0.5))
            iy = int(np.floor(y + 0.5))
            ix_arr[i] = ix
            iy_arr[i] = iy
            offset = galsim.PositionD(x - ix, y - iy)

            # Create the sub-image for this galaxy
            if config.get("skymodel", "galaxy_profile") == "points":
                stamp = gal.drawImage(
                    scale=pixel_scale / galsim.arcsec,
                    offset=offset,
                    method="real_space",
                )
            else:
                stamp = gal.drawImage(scale=pixel_scale / galsim.arcsec, offset=offset)

            stamp.setCenter(ix, iy)

            # Add the sub-image to the full iamge
            bounds = stamp.bounds & full_image.bounds
            full_image[bounds] += stamp[bounds]
            sys.stdout.flush()

    if config.getboolean("skymodel", "doagn"):
        # Load the catalogue
        cat_agn = Table()
        cat_file_name = config.get("skymodel", "agn_catalogue_filepath")
        print("Loading catalogue from {0} ...".format(cat_file_name))
        cat_read_agn = Table.read(cat_file_name, format="ascii")

        source_prefix = "TRECS-"
        source_name_ra = np.asarray(cat_read_agn["lon(deg)"], dtype=str)
        source_name_dec = np.asarray(cat_read_agn["lat(deg)"], dtype=str)

        source_prefix_arr = strmultiply(
            source_prefix, np.ones_like(source_name_ra, dtype=int)
        )
        source_l_arr = strmultiply("l", np.ones_like(source_name_ra, dtype=int))
        source_b_arr = strmultiply("b", np.ones_like(source_name_ra, dtype=int))

        source_name_pos = stradd(
            source_prefix_arr,
            stradd(
                source_l_arr,
                stradd(source_name_ra, (stradd(source_b_arr, source_name_dec))),
            ),
        )
        cat_agn["Source_id"] = source_name_pos

        cat_agn["ra_offset"] = cat_read_agn["lon(deg)"]  # deg
        cat_agn["ra_offset"].unit = "deg"

        cat_agn["dec_offset"] = cat_read_agn["lat(deg)"]  # deg
        cat_agn["dec_offset"].unit = "deg"

        cat_agn["DEC"] = dec_field_gs / galsim.degrees + cat_agn["dec_offset"]
        dec_abs_radians = cat_agn["DEC"] * galsim.degrees / galsim.radians
        cat_agn["RA"] = ra_field_gs / galsim.degrees + cat_agn["ra_offset"] / np.cos(
            np.asarray(dec_abs_radians, dtype=float)
        )

        cat_agn["Total_flux"] = cat_read_agn["flux(mJy)"] * 1.0e-3  # Jy
        cat_agn["Total_flux"].unit = "Jy"

        cat_agn["Maj"] = cat_read_agn["size(arcsec)"]  # arcsec
        cat_agn["Maj"].unit = "arcsec"

        cat_agn["e1"] = cat_read_agn["e1"]
        cat_agn["e2"] = cat_read_agn["e2"]

        cat_agn["g1_shear"] = cat_read_agn["gamma1"]
        cat_agn["g2_shear"] = cat_read_agn["gamma2"]

        cat_agn["mod_e"] = np.sqrt(cat_agn["e1"] ** 2.0 + cat_agn["e2"] ** 2.0)

        cat_agn["Rs"] = cat_read_agn["Rs"]

        cat_agn["PA"] = random.uniform(0, 2.0 * np.pi, len(cat_read_agn))
        cat_agn["PA"].unit = "rad"

        # fov cut
        ra_offset_max = 0.9 * (fov / 2) / galsim.degrees
        dec_offset_max = 0.9 * (fov / 2) / galsim.degrees

        fov_cut = (abs(cat_agn["ra_offset"]) < ra_offset_max) * (
            abs(cat_agn["dec_offset"]) < dec_offset_max
        )
        cat_agn = cat_agn[fov_cut]

        # flux cuts
        if config.getboolean("skymodel", "highfluxcut"):
            highflux_cut = cat_agn["Total_flux"] < config.getfloat(
                "skymodel", "highfluxcut_value"
            )
            cat_agn = cat_agn[highflux_cut]

        if config.getboolean("skymodel", "lowfluxcut"):
            lowflux_cut = cat_agn["Total_flux"] > config.getfloat(
                "skymodel", "lowfluxcut_value"
            )
            cat_agn = cat_agn[lowflux_cut]

        if config.getboolean("skymodel", "highsizecut"):
            highsize_cut = cat_agn["Maj"] < config.getfloat(
                "skymodel", "highsizecut_value"
            )
            cat_agn = cat_agn[highsize_cut]

        if config.getboolean("skymodel", "lowsizecut"):
            lowsize_cut = cat_agn["Maj"] > config.getfloat(
                "skymodel", "lowsizecut_value"
            )
            cat_agn = cat_agn[lowsize_cut]

        # number of sources, on grid if requested
        if config.getboolean("skymodel", "grid"):
            nobj = int(np.sqrt(config.getint("skymodel", "ngals"))) ** 2.0
            cat_agn["ra_offset"] = np.linspace(-ra_offset_max, ra_offset_max, nobj)
            cat_agn["dec_offset"] = np.linspace(-ra_offset_max, ra_offset_max, nobj)
        else:
            nobj = len(cat_agn)
            if config.getint("skymodel", "ngals") > -1:
                nobj = config.getint("skymodel", "ngals")
                cat_agn = cat_agn[:nobj]

        # flux range
        if config.get("skymodel", "fluxscale") == "constant":
            cat_agn["Total_flux"] = np.ones_like(
                cat_agn["Total_flux"]
            ) * config.getfloat("skymodel", "fluxscale_constant_value")

        # scale flux
        cat_agn["Total_flux"] = cat_agn["Total_flux"] * config.getfloat(
            "skymodel", "flux_factor"
        )

        # scale size
        cat_agn["Maj"] = cat_agn["Maj"] * config.getfloat("skymodel", "sizefactor")

        # write out catalogue

        output_filename = (
            "agn_truthcat" + config.get("pipeline", "output_suffix") + ".fits"
        )

        print("Writing truthcat to " + os.path.join(output_path, output_filename))
        cat.write(
            os.path.join(output_path, output_filename), format="fits", overwrite=True,
        )

        ix_arr = np.ones(nobj)
        iy_arr = np.ones(nobj)
        print("...done.")

        # Draw the galaxies onto the galsim image
        for i, cat_gal in enumerate(cat_agn):

            sys.stdout.write(
                "\rAdding agn source {0} of {1} to skymodel...".format(i + 1, nobj)
            )

            if (cat_gal["Rs"] < 0.01) or (
                cat_gal["Maj"] < config.getfloat("skymodel", "pixel_scale") / 2
            ):
                x, y = w_twod.wcs_world2pix(cat_gal["RA"], cat_gal["DEC"], 0,)
                x = float(x)
                y = float(y)

                # Account for the fractional part of the position:
                ix = int(np.floor(x + 0.5))
                iy = int(np.floor(y + 0.5))
                ix_arr[i] = ix
                iy_arr[i] = iy
                offset = galsim.PositionD(x - ix, y - iy)

                # Create the sub-image for this galaxy
                gal = galsim.Gaussian(fwhm=1.0, flux=1.0e-10)
                stamp = gal.drawImage(scale=pixel_scale / galsim.arcsec, offset=offset)
                stamp.setCenter(ix, iy)

                cen = stamp.array.shape
                # Add the hotspots as single pixel point sources
                stamp.array[cen[0] / 2, cen[1] / 2] += cat_gal["Total_flux"]

                # Add the sub-image to the full iamge
                bounds = stamp.bounds & full_image.bounds
                full_image[bounds] += stamp[bounds]

            else:
                lobe_flux = cat_gal["Total_flux"] * 0.99
                hs_flux = cat_gal["Total_flux"] - lobe_flux
                hs1_flux = hs_flux / 3.0
                hs2_flux = hs_flux / 3.0
                hs3_flux = hs_flux / 3.0

                hs_offset = cat_gal["Rs"] * cat_gal["Maj"]
                lobe_offset = cat_gal["Maj"] * 0.6

                lobe1 = galsim.Gaussian(
                    sigma=cat_gal["Maj"] * 0.25,
                    flux=lobe_flux / 2.0,
                    gsparams=big_fft_params,
                )
                lobe2 = galsim.Gaussian(
                    sigma=cat_gal["Maj"] * 0.25,
                    flux=lobe_flux / 2.0,
                    gsparams=big_fft_params,
                )

                lobe1 = lobe1.shear(e1=0.3, e2=0)
                lobe2 = lobe2.shear(e1=0.3, e2=0)

                lobe1 = lobe1.shift(-lobe_offset, 0)
                lobe2 = lobe2.shift(lobe_offset, 0)

                gal = lobe1 + lobe2
                gal = gal.rotate(cat_gal["PA"] * galsim.radians)

                x, y = w_twod.wcs_world2pix(cat_gal["RA"], cat_gal["DEC"], 0,)
                x = float(x)
                y = float(y)

                # Account for the fractional part of the position:
                ix = int(np.floor(x + 0.5))
                iy = int(np.floor(y + 0.5))
                ix_arr[i] = ix
                iy_arr[i] = iy
                offset = galsim.PositionD(x - ix, y - iy)
                hs_offset_pixels = hs_offset * pixel_scale / galsim.arcsec
                hs_ix_offset = hs_offset * np.sin(cat_gal["PA"])
                hs_iy_offset = hs_offset * np.cos(cat_gal["PA"])
                hs_ix_offset_pixels = int(hs_ix_offset * pixel_scale / galsim.arcsec)
                hs_iy_offset_pixels = int(hs_iy_offset * pixel_scale / galsim.arcsec)

                # Create the sub-image for this galaxy
                stamp = gal.drawImage(scale=pixel_scale / galsim.arcsec, offset=offset)
                stamp.setCenter(ix, iy)

                cen = stamp.array.shape
                # Add the hotspots as single pixel point sources
                stamp.array[cen[0] / 2, cen[1] / 2] += hs1_flux

                stamp.array[
                    cen[0] / 2 + hs_ix_offset_pixels, cen[1] / 2 + hs_iy_offset_pixels
                ] += hs2_flux
                stamp.array[
                    cen[0] / 2 - hs_ix_offset_pixels, cen[1] / 2 - hs_iy_offset_pixels
                ] += hs3_flux

                if config.getboolean("skymodel", "pickleagn"):
                    pickle.dump(stamp.array, open("agn_{0}.p".format(i), "wb"))

                # Add the sub-image to the full iamge
                bounds = stamp.bounds & full_image.bounds
                full_image[bounds] += stamp[bounds]

            sys.stdout.flush()

    tend = time.time()
    print("\n...done in {0} seconds.".format(tend - tstart))

    output_image_filename = config.get("pipeline", "output_suffix") + "_skymodel.fits"

    print(
        "Writing {1} image data to {0} ...".format(
            os.path.join(output_path, output_image_filename), full_image.array.shape
        )
    )

    config.set(
        "observation",
        "image_filepath",
        os.path.join(output_path, output_image_filename),
    )

    # Extract the numpy array from the galsim image
    image_data = full_image.array

    hdulist = fits.HDUList([hdu])
    hdulist.writeto(os.path.join(output_path, output_image_filename), clobber=True)

    print("...done.")

    if config.getboolean("skymodel", "make_im3cat"):
        output_im3cat_filename = config.get("pipeline", "output_suffix") + "_im3cat.txt"

        np.savetxt(
            os.path.join(output_path, output_im3cat_filename),
            np.column_stack([np.arange(nobj), ix_arr, iy_arr]),
        )

    print("runSkyModel complete.")
