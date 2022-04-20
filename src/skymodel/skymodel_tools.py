import numpy as np
import pdb
from astropy.io import fits
from astropy import wcs
import galsim


def setup_wcs(config, ndim, nu_axis=False):

    pixel_scale = config.getfloat("skymodel", "pixel_scale") * galsim.arcsec
    fov = config.getfloat("skymodel", "field_of_view") * galsim.arcmin
    image_size = int((fov / galsim.arcmin) / (pixel_scale / galsim.arcmin))

    ra_field = config.get("observation", "field_ra")
    ra_field_gs = galsim.HMS_Angle(ra_field)
    dec_field = config.get("observation", "field_dec")
    dec_field_gs = galsim.DMS_Angle(dec_field)

    n_ifs = config.getint("observation", "n_IFs")
    bw = config.getfloat("observation", "total_bandwidth")
    base_freq = config.getfloat("observation", "lowest_frequency")
    n_chan = config.getint("observation", "n_channels")

    channel_width = bw / (n_chan * n_ifs)
    if_width = bw / n_ifs

    w = wcs.WCS(naxis=ndim)
    if ndim == 4:
        if nu_axis:
            w.wcs.crpix = [float(image_size) / 2, float(image_size) / 2, 1, 1]
            w.wcs.cdelt = [
                pixel_scale / galsim.degrees,
                pixel_scale / galsim.degrees,
                bw,
                1,
            ]
            w.wcs.crval = [
                ra_field_gs / galsim.degrees,
                dec_field_gs / galsim.degrees,
                base_freq + bw / 2,
                1,
            ]
            w.wcs.ctype = ["RA---SIN", "DEC--SIN", "FREQ", "STOKES"]
            w.wcs.cunit = ["deg", "deg", "Hz", ""]

        else:
            w.wcs.crpix = [float(image_size) / 2, float(image_size) / 2, 1, 1]
            w.wcs.cdelt = [
                pixel_scale / galsim.degrees,
                pixel_scale / galsim.degrees,
                bw,
                1,
            ]
            w.wcs.crval = [
                ra_field_gs / galsim.degrees,
                dec_field_gs / galsim.degrees,
                base_freq + bw / 2,
                1,
            ]
            w.wcs.ctype = ["RA---SIN", "DEC--SIN", "FREQ", "STOKES"]
            w.wcs.cunit = ["deg", "deg", "Hz", ""]

    elif ndim == 2:
        w.wcs.crpix = [float(image_size) / 2, float(image_size) / 2]
        w.wcs.cdelt = [pixel_scale / galsim.degrees, pixel_scale / galsim.degrees]
        w.wcs.crval = [ra_field_gs / galsim.degrees, dec_field_gs / galsim.degrees]
        w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
        w.wcs.cunit = ["deg", "deg"]

    return w


def write4dImage(
    outname,
    image_data,
    pixel_scale,
    obs_ra=150.0e0,
    obs_dec=90.0e0,
    obs_freq=1.4e9,
    obs_bwidth=1.0e9,
    clobber=True,
):
    """Write a FITS file with minimum necessary information to use as input for a
  CASA observation simulation.

  Parameters
  ----------
  outname : string
    Output file name and location.
  image_data : ndarray
    Numpy array containing image data. Currently expects a 2d array.
  pixel_scale : float
    Scale of pixels in image ***IN DEGREES***
  obs_ra : float, optional, default: 150.e0
    Right Ascension of field **IN DEGREES***
  obs_dec : float, optional, default: 90.e0
    Declination of field ***IN DEGREES***
  obs_freq : float, optional, default: 1.4e9
    Frequency of emission in image ***in Hz***
  obs_bwidth : float, optional, default: 1.e9
    Bandwidth of image ***in Hz***. Should cover entire expected
    bandwidth of simulated observation.

  """
    image_dim = image_data.ndim
    # TODO: some error checking on ndim (i.e. if already 4)

    new_image = np.zeros([1, 1, image_data.shape[0], image_data.shape[1]])
    new_image[0, 0, :, :] = image_data
    new_image = np.asarray(new_image, dtype=np.float32)

    newhdr = fits.Header()

    newhdr["SIMPLE"] = True
    newhdr["BITPIX"] = -32
    newhdr["NAXIS"] = 4
    newhdr["NAXIS1"] = image_data.shape[0]
    newhdr["NAXIS2"] = image_data.shape[1]
    newhdr["NAXIS3"] = 1
    newhdr["NAXIS4"] = 1

    newhdr["EXTEND"] = True
    newhdr["BSCALE"] = 1.0e0
    newhdr["BZERO"] = 0.0e0
    newhdr["BMAJ"] = pixel_scale
    newhdr["BMIN"] = pixel_scale
    newhdr["BTYPE"] = "Intensity"
    newhdr["OBJECT"] = "        "

    newhdr["BUNIT"] = "JY/PIXEL"
    newhdr["EQUINOX"] = 2.0e3
    newhdr["RADESYS"] = "FK5"
    newhdr["LONPOLE"] = 1.8e2
    newhdr["LATPOLE"] = 6.8e1
    newhdr["PC01_01"] = 1.000000000000e00
    newhdr["PC02_01"] = 0.000000000000e00
    newhdr["PC03_01"] = 0.000000000000e00
    newhdr["PC04_01"] = 0.000000000000e00
    newhdr["PC01_02"] = 0.000000000000e00
    newhdr["PC02_02"] = 1.000000000000e00
    newhdr["PC03_02"] = 0.000000000000e00
    newhdr["PC04_02"] = 0.000000000000e00
    newhdr["PC01_03"] = 0.000000000000e00
    newhdr["PC02_03"] = 0.000000000000e00
    newhdr["PC03_03"] = 1.000000000000e00
    newhdr["PC04_03"] = 0.000000000000e00
    newhdr["PC01_04"] = 0.000000000000e00
    newhdr["PC02_04"] = 0.000000000000e00
    newhdr["PC03_04"] = 0.000000000000e00
    newhdr["PC04_04"] = 1.000000000000e00
    newhdr["CTYPE1"] = "RA---SIN"
    newhdr["CRVAL1"] = obs_ra
    newhdr["CDELT1"] = pixel_scale
    newhdr["CRPIX1"] = float(new_image.shape[2] / 2)
    newhdr["CUNIT1"] = "deg"

    newhdr["CTYPE2"] = "DEC--SIN"
    newhdr["CRVAL2"] = obs_dec
    newhdr["CDELT2"] = pixel_scale
    newhdr["CRPIX2"] = float(new_image.shape[3] / 2)
    newhdr["CUNIT2"] = "deg"

    newhdr["CTYPE3"] = "STOKES"
    newhdr["CRVAL3"] = 1.0e0
    newhdr["CDELT3"] = 1.0e0
    newhdr["CRPIX3"] = 1.0e0
    newhdr["CUNIT3"] = "        "

    newhdr["CTYPE4"] = "FREQ"
    newhdr["CRVAL4"] = obs_freq
    newhdr["CDELT4"] = obs_bwidth
    newhdr["CRPIX4"] = 1.0e0
    newhdr["CUNIT4"] = "Hz      "

    newhdr["PV2_1"] = 0.0e0
    newhdr["PV2_2"] = 0.0e0

    newhdr["RESTFRQ"] = obs_freq
    newhdr["SPECSYS"] = "LSRK    "
    newhdr["ALTRVAL"] = -0.0e0
    newhdr["ALTRPIX"] = 1.000000000000e00
    newhdr["VELREF"] = 257
    newhdr["TELESCOP"] = "SKYMODEL "
    newhdr["OBSERVER"] = "A. R. Bitrary"
    newhdr["DATE-OBS"] = "2000-01-01T00:00:00.000100"
    newhdr["TIMESYS"] = "UTC     "
    newhdr["OBSRA"] = obs_ra
    newhdr["OBSDEC"] = obs_dec
    newhdr["OBSGEO-X"] = 2.225061873184e06
    newhdr["OBSGEO-Y"] = -5.440061952280e06
    newhdr["OBSGEO-Z"] = -2.481682085791e06
    newhdr["DATE"] = "2014-02-20T13:17:20.141840"
    newhdr["ORIGIN"] = "write4dImage"

    hdu = fits.PrimaryHDU(new_image, header=newhdr)
    hdu.header.insert(8, ("EXTEND", "BSCALE", "BZERO"))
    hdu.header["EXTEND"] = True
    hdu.header["BSCALE"] = 1.0e0
    hdu.header["BZERO"] = 0.0e0
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(outname, clobber=clobber)


def ellip_pdf(x, B=0.19, C=0.58):
    """Ellipticity probability distribution function from GREAT08
  (arXiv:0908.0945 equation A3).

  Parameters
  ----------
  x : float
    Ellipticity modulus at which to calculate PDF
  B : float
    Scale parameter
  C : float
    Index parameter
  """
    retVar = x * (np.cos(np.pi * x / 2.0) ** 2.0) * np.exp(-((2.0 * x / B) ** C))
    return retVar
