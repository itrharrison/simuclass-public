import ConfigParser
import cPickle as pickle
import pdb
import os
from astropy.io import fits


def runCASAClean(config, parameter_filename):

    run_dir = config.get("pipeline", "simuclass_dirname")

    os.chdir(config.get("pipeline", "output_dirname"))
    temp_dir = os.getcwd()
    pickle.dump(config, open(temp_dir + "/temp_config.p", "wb"))

    casa_exe = config.get("pipeline", "casa_exe")
    cmd = "{3} --nogui --log2term -c {0}/src/imager/casa_imager.py {1}/{2}".format( \
        run_dir, temp_dir, "temp_config.p", casa_exe
    )
    os.system(cmd)

    os.chdir(run_dir)

def runWsclean(config, parameter_filename):

    run_dir = os.getcwd()

    os.chdir(config.get("pipeline", "output_dirname"))

    ms_list = ' '.join(config.get('observation', 'uvcoverage_filename').split(','))

    wsclean_exe = config.get("pipeline", "wsclean_exe")
    
    cmd = 'wsclean -j ' + config.get('imager', 'ncores')\
            '-abs-mem ' + config.get('imager', 'memory_size')\
            '-grid-with-beam'\
            '-use-idg'\
            '-idg-mode cpu'\
            '-aterm-kernel-size ' + config.get('imager', 'aterm_kernel_size')\
            '-name ' + config.get('pipeline', 'output_suffix')\
            '-size ' + config.get('imager', 'npix') + ' ' + config.get('imager', 'npix')\
            '-scale ' + + config.get('imager', 'pixel_scale') + 'asec'\
            '-mgain ' + config.get('imager', 'm_gain')\
            '-weight briggs ' + config.get('imager', 'briggs_robust')\
            '-data-column DATA'\
            '-auto-mask ' + config.get('imager', 'auto_masking')\
            '-auto-threshold ' + config.get('imager', 'auto_threshold')\
            '-niter ' + config.get('imager', 'number_of_iterations')\
            '-parallel-deconvolution ' + config.get('imager', 'parallel_size')\
            '-no-update-model-required'\
            '-minuv-l ' + config.get('imager', 'min_uv_l')\
            '-maxuv-l ' + config.get('imager', 'max_uv_l')\
            '-taper-inner-tukey ' + config.get('imager', 'tukey_inner_taper')\
            '-taper-tukey ' + config.get('imager', 'tukey_taper')\
            '-log-time'\
            + ms_list

    os.system(cmd)

    os.chdir(run_dir)

