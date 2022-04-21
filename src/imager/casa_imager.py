import ConfigParser
import cPickle as pickle
import pdb
import os
import numpy as np

if __name__ == "__main__":

    config = pickle.load(open(sys.argv[-1], "rb"))

    outdir = config.get("pipeline", "output_dirname")
    msdir = config.get("observation", "uvcoverage_dirname")
    msname = config.get("observation", "uvcoverage_filename")
    suffix = config.get("pipeline", "output_suffix")

    imagename = msname + "_" + suffix + ".tclean"
    npix = config.getfloat("imager", "npix")
    pixel_scale = config.getfloat("imager", "pixel_scale")
    niter_value = config.getint("imager", "number_of_iterations")
    threshold = config.getfloat("imager", "threshold")
    deconvolver = config.get("imager", "deconvolver")
    beam_maj = config.get("imager", "beam_maj")
    beam_min = config.get("imager", "beam_min")
    beam_pa = config.get("imager", "beam_pa")
    pbcor = config.getboolean("imager", "pbcor")
    scales = config.get("imager", "multi_scales").split(",")
    scales = np.asarray(scales, dtype=int)
    nterms = config.getint("imager", "nterms")
    weighting = config.get("imager", "weighting")
    robust = config.get("imager", "briggs_robust")
    cell_size = [str(pixel_scale) + "arcsec"]
    reffreq = config.get("imager", "ref_freq")
    export_only = config.getboolean("imager", "export_only")
    restbeam = [beam_maj + "arcsec", beam_min + "arcsec", beam_pa + "deg"]

    if not export_only:
        tclean(
            vis=os.path.join(msdir, msname),
            imagename=os.path.join(outdir, imagename),
            threshold=threshold,
            niter=niter_value,
            imsize=[int(npix), int(npix)],
            cell=cell_size,
            deconvolver=deconvolver,
            scales=scales,
            nterms=nterms,
            weighting=weighting,
            pbcor=pbcor,
            pblimit=0.00001,
            reffreq=reffreq,
        )

    exportfits(
        imagename=os.path.join(outdir, imagename + ".image.tt0"),
        fitsimage=os.path.join(outdir, imagename + ".image.tt0.fits"),
        overwrite=True,
    )

    exportfits(
        imagename=os.path.join(outdir, imagename + ".residual.tt0"),
        fitsimage=os.path.join(outdir, imagename + ".residual.tt0.fits"),
        overwrite=True,
    )
