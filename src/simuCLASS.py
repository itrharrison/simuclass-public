import pdb
import sys
import os
import ConfigParser

config = ConfigParser.ConfigParser()
config.read(sys.argv[1])

if not os.path.exists(config.get("pipeline", "output_dirname")):
    os.mkdir(config.get("pipeline", "output_dirname"))

if config.getboolean("pipeline", "doskymodel"):
    from skymodel.skymodel import runSkyModel

    print("Running skymodel...")
    runSkyModel(config)

if config.getboolean("pipeline", "dosimdata"):
    from observation.observation import runSimulateData

    print("Running simulatedata...")
    runSimulateData(config)

if config.getboolean("pipeline", "doimagedata"):

    if config.get("imager", "imager_type") == 'wsclean':

        from imager.imager import runWsclean as runImager
        print("Running wsclean imager...")

    else:

        from imager.imager import runCASAClean as runImager
        print("Running CASA tclean imager...")

    runImager(config, sys.argv[1])
