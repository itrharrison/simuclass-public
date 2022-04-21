import numpy as np
import cPickle as pickle
import os
import time
import pdb
import ConfigParser


def runSimulateData(config):

    run_dir = os.getcwd()

    os.chdir(config.get("pipeline", "output_dirname"))
    temp_dir = os.getcwd()
    config.set("pipeline", "run_dir", os.path.dirname(os.path.abspath(__file__)))
    pickle.dump(config, open(temp_dir + "/temp_config.p", "wb"))

    casa_exe = config.get("pipeline", "casa_exe")
    cmd = "{3} --nogui --log2term -c {0}/src/observation/observation.py {1}/{2}".format(
        run_dir, temp_dir, "temp_config.p", casa_exe
    )
    print (cmd)
    os.system(cmd)

    os.chdir(run_dir)


if __name__ == "__main__":


    config = pickle.load(open(sys.argv[-1], "rb"))

    sys.path.append(config.get('pipeline', 'run_dir'))
    from emerlin_antenna import *
    from vla_antenna import *

    n_ifs = config.getint("observation", "n_IFs")
    bw = config.getfloat("observation", "total_bandwidth")
    base_freq = config.getfloat("observation", "lowest_frequency")
    channel_width = config.getfloat("observation", "channel_width")
    n_chan = config.getint("observation", "n_channels")
    channel_separation = config.getfloat("observation", "channel_separation")
    output_path = config.get("pipeline", "output_dirname")

    if_width = bw / n_ifs

    if config.get("observation", "uvcoverage_type") == "simulate":
        msname = config.get("pipeline", "output_suffix") + ".ms"
        msname = os.path.join(output_path, msname)
        print ("making uv coverage! \n at {0}".format(msname))
        sm.open(msname)
        for i in np.arange(1, n_ifs + 1):
            fr = base_freq + if_width * (i - 1)
            print (str(channel_width) + "Hz")
            sm.setspwindow(
                spwname="IF" + str(i),
                freq=str(fr) + "Hz",  # starting frequency
                deltafreq=str(channel_separation) + "Hz",  # increment per chan
                freqresolution=str(channel_width) + "Hz",  # width per chan
                nchannels=n_chan,
                stokes="I",
            )

        if config.get("observation", "telescope") == "e-merlin":
            observatory = "e-MERLIN"
            posemerlin = me.observatory(observatory)
            sm.setconfig(
                telescopename="e-MERLIN",
                x=emerlin_xx,
                y=emerlin_yy,
                z=emerlin_zz,
                dishdiameter=emerlin_diam.tolist(),
                mount="alt-az",
                coordsystem="global",
            )
            sm.setfeed(mode="perfect R L")
        elif config.get("observation", "telescope") == "jvla":
            posvla = me.observatory("vla")
            sm.setconfig(
                telescopename="VLA",
                x=vla_xx,
                y=vla_yy,
                z=vla_zz,
                dishdiameter=vla_diam.tolist(),
                mount="alt-az",
                coordsystem="global",
            )
            sm.setfeed(mode="perfect R L")

        source_dec_casa = config.get("observation", "field_dec").split(":")
        source_dec_casa = (
            source_dec_casa[0]
            + "d"
            + source_dec_casa[1]
            + "m"
            + source_dec_casa[2]
            + "s"
        )
        source_dirn = me.direction(
            "J2000", config.get("observation", "field_ra"), source_dec_casa
        )
        sm.setfield(
            sourcename=config.get("observation", "field_name"),
            sourcedirection=source_dirn,
        )

        obs_date = time.strftime("%Y/%m/%d", time.gmtime())
        ref_time = me.epoch("IAT", obs_date)
        sm.settimes(
            integrationtime=config.get("observation", "t_int") + "s",
            usehourangle=True,
            referencetime=ref_time,
        )

        for i in np.arange(1, n_ifs + 1):
            sm.observe(
                config.get("observation", "field_name"),
                "IF" + str(i),
                starttime="0s",
                stoptime=config.get("observation", "observation_time") + "s",
            )
    else:
        print ("loading uv coverage!")
        msdir = config.get("observation", "uvcoverage_dirname")
        msname = config.get("observation", "uvcoverage_filename")
        sm.openfromms(os.path.join(msdir, msname))

    fitsimage = config.get("observation", "image_filepath")
    imagename = fitsimage.replace(".fits", ".im")

    importfits(fitsimage=fitsimage, imagename=imagename, overwrite=True)

    sm.predict(imagename=imagename)

    if config.get("observation", "noisemode") == "uniform":
        sm.setnoise(
            mode="simplenoise",
            simplenoise=config.get("observation", "uniform_noise") + "Jy",
        )
        sm.corrupt()
    elif config.get("observation", "noisemode") == "none":
        pass
    elif config.get("observation", "noisemode") == "real":
        rms_noise_list = pickle.load(
            open(config.get("observation", "noise_filepath"), "rb")
        )
        for bl_spw in rms_noise_list:
            sm.openfromms(msname)
            # this_selection = {'baseline' : bl_spw[0], 'spw': bl_spw[1]}
            print bl_spw
            ant1 = bl_spw[0].split("&")[0]
            ant2 = bl_spw[0].split("&")[1]

            # CASA antenna naming schemes are horrible and inconsistent!
            # Convert from Antenna 'Name' to 'ID'
            ant1 = ant1.replace("1", "0")
            ant1 = ant1.replace("2", "1")
            ant1 = ant1.replace("5", "4")
            ant1 = ant1.replace("6", "5")
            ant1 = ant1.replace("7", "6")
            ant1 = ant1.replace("8", "7")
            ant1 = ant1.replace("9", "8")

            ant2 = ant2.replace("1", "0")
            ant2 = ant2.replace("2", "1")
            ant2 = ant2.replace("5", "4")
            ant2 = ant2.replace("6", "5")
            ant2 = ant2.replace("7", "6")
            ant2 = ant2.replace("8", "7")
            ant2 = ant2.replace("9", "8")

            this_selection = "ANTENNA1=={0} and ANTENNA2=={1}".format(ant1, ant2)
            this_rms = rms_noise_list[bl_spw]
            print this_selection
            print bl_spw[1]
            print str(this_rms) + "Jy"
            sm.setdata(spwid=int(bl_spw[1]), msselect=this_selection)
            sm.setnoise(mode="simplenoise", simplenoise=str(this_rms) + "Jy")
            sm.corrupt()
            sm.done()
