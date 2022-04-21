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
    cmd = "{3} --nogui --log2term -c {0}/src/imager/casa_imager.py {1}/{2}".format(
        run_dir, temp_dir, "temp_config.p", casa_exe
    )
    os.system(cmd)

    os.chdir(run_dir)
