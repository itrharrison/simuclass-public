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

    runSkyModel(config)

if config.getboolean("pipeline", "dosimdata"):
    from observation.observation import runSimulateData

    runSimulateData(config)

if config.getboolean("pipeline", "doimagedata"):
    from imager.imager import runCASAClean

    runCASAClean(config, sys.argv[1])
