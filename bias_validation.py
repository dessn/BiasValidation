#!/usr/bin/env python
import multiprocessing
from multiprocessing.pool import Pool
import argparse
import logging
import coloredlogs
import tomli as tml
from tqdm import tqdm
from pathlib import Path
import difflib
import shutil
import subprocess
import sys
import os
from files.scripts.analysis import *

config = {
    "logfile": "log.txt",
    "base": "/project2/rkessler/SURVEYS/DES/USERS/parmstrong/dev/bias_validation/",
    "clean": True
}

def setup_logging(logging_filename, verbose):
    level = logging.DEBUG if verbose else logging.INFO

    fmt = "[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s" if verbose else "%(message)s"
    handlers = [logging.StreamHandler()]
    handlers[0].setLevel(level)
    handlers.append(logging.FileHandler(logging_filename, mode="w"))
    handlers[-1].setLevel(logging.DEBUG)
    logging.basicConfig(level=level, format=fmt, handlers=handlers)

    coloredlogs.install(
        level = level,
        fmt = fmt,
        reconfigure = True,
        level_styles = coloredlogs.parse_encoded_styles("debug=8;info=green;warning=yellow;error=red,bold;critical=red,inverse"),
    )
    logging.getLogger("matplotlib").setLevel(logging.ERROR)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("toml", help="the name of the toml config file to run.", type=str)
    parser.add_argument("-v", "--verbose", help="Increase output verbosity", action="store_true")
    parser.add_argument("-n", "--no_skip", help="Rerun pippin jobs", action="store_true")

    args = parser.parse_args()
    return args

def validation_setup(options, config, logger):
    if options is None:
        logger.error("No validation cosmology specified, quitting")
        return None
    name = list(options.keys())[0]
    options = options[name]
    Om_list = options["OMEGA_MATTER"]
    w0_list = options["W0_LAMBDA"]
    input_path = config.get("pippin_inputs") / f"BV_V_{name}.yml"
    master = config.get("input_files") / "validation_master.yml"
    # Edit the master yml file
    SIM_MASTER = """DESSIM_{i}:
        IA_G10_DES3YR:
            BASE: /project2/rkessler/SURVEYS/DES/USERS/parmstrong/dev/bias_validation/files/base_files/sn_ia_salt2_g10_des3yr.input
        GLOBAL:
            W0_LAMBDA: {W0_LAMBDA}
            OMEGA_MATTER: {OMEGA_MATTER}
            OMEGA_LAMBDA: {OMEGA_LAMBDA}
            NGEN_UNIT: 5.0
            RANSEED_CHANGE: 150 12345
    DESBIASCOR_{i}:
        IA_G10_DES3YR:
            BASE: /project2/rkessler/SURVEYS/DES/USERS/parmstrong/dev/bias_validation/files/base_files/sn_ia_salt2_g10_des3yr.input
            GENSIGMA_SALT2ALPHA: 1E8  1E8
            GENRANGE_SALT2ALPHA: 0.12 0.20
            GENGRID_SALT2ALPHA: 2
            GENSIGMA_SALT2BETA: 1E8  1E8
            GENRANGE_SALT2BETA: 2.6  3.6
            GENGRID_SALT2BETA: 2
        GLOBAL:
            W0_LAMBDA: {W0_LAMBDA}
            OMEGA_MATTER: {OMEGA_MATTER}
            OMEGA_LAMBDA: {OMEGA_LAMBDA}
            NGEN_UNIT: 50.0
            RANSEED_REPEAT: 50 12345
    """
    REPLACE_SIM = ""
    for i in range(len(Om_list)):
        Om = Om_list[i]
        w0 = w0_list[i]
        SIM = SIM_MASTER.format(**{"i": i, "OMEGA_MATTER": Om, "OMEGA_LAMBDA": 1-Om, "W0_LAMBDA": w0})
        REPLACE_SIM += SIM
    BIASCOR_MASTER = """BCOR_OMW_NO_OMPRI_{i}:
        BASE: /project2/rkessler/SURVEYS/DES/USERS/parmstrong/dev/bias_validation/files/base_files/SALT2mu_no_ompri.input
        DATA: [DESFITSYS_DESSIM_{i}]
        SIMFILE_BIASCOR: [DESFIT_DESBIASCOR_{i}]
        NMAX: 190
    """
    REPLACE_BIASCOR = ""
    for i in range(len(Om_list)):
        Om = Om_list[i]
        w0 = w0_list[i]
        BIASCOR = BIASCOR_MASTER.format(**{"i": i, "OMEGA_MATTER": Om, "W0_LAMBDA": w0})
        REPLACE_BIASCOR += BIASCOR 
    COV_MASTER = """COV_OMW_NO_OMPRI_{i}:
        MASK: BCOR_OMW_NO_OMPRI_{i}
    """
    REPLACE_COV = ""
    for i in range(len(Om_list)):
        COV = COV_MASTER.format(**{"i": i})
        REPLACE_COV += COV
    WFIT_MASTER = """   SN_OMW_NO_OMPRI_{i}:
            MASK: COV_OMW_NO_OMPRI_{i}
            OPTS:
                BATCH_INFO: sbatch $SBATCH_TEMPLATES/SBATCH_Midway2b.TEMPLATE 50
                WFITOPTS:
                    - /om_no_pri/ -ompri {Om} -dompri 50 -outfile_chi2grid chi2.txt -minchi2 -omsteps 201 -ommin {Om_min} -ommax {Om_max} -w0steps 301 -w0min {w0_min} -w0max {w0_max}
                WFITAVG:
                    - COV_OMW_NO_OMPRI_{i}_BCOR_OMW_NO_OMPRI_{i}
    """
    REPLACE_WFIT = ""
    for i in range(len(Om_list)):
        Om = Om_list[i]
        Om_min = Om - 1
        Om_max = Om + 1
        w0 = w0_list[i]
        w0_min = w0 - 2
        w0_max = w0 + 2
        WFIT = WFIT_MASTER.format(**{"i": i, "Om": Om, "w0": w0, "Om_min": Om_min, "Om_max": Om_max, "w0_min": w0_min, "w0_max": w0_max})
        REPLACE_WFIT += WFIT
    with open(master, 'r') as f:
        tmp = f.read()
    prepend = f"# Generated via:\n# [ validation.{name} ]\n# OMEGA_MATTER = {Om_list}\n# W0_LAMBDA = {w0_list}"
    tmp = tmp.format(**{"PREPEND": prepend, "REPLACE_SIM": REPLACE_SIM, "REPLACE_BIASCOR": REPLACE_BIASCOR, "REPLACE_COV": REPLACE_COV, "REPLACE_WFIT": REPLACE_WFIT})
    skip = False
    if input_path.exists():
        logger.info(f"Found another input file under {name}, comparing them now")
        with open(input_path, 'r') as f:
            tmp2 = f.read()
        diffs = list(difflib.unified_diff(tmp.split('\n'), tmp2.split('\n'), fromfile="new_input", tofile="old_input", lineterm=''))
        if len(diffs) == 0:
            logger.info("No differences found, will not rerun")
            skip = True
        else:
            logger.warning(f"Found the following differences:")
            for line in diffs:
                logger.warning(line)
            logger.warning("If you continue, you will overwrite previously defined inputs. Would you like to continue (y/n)?")
            response = input().upper()
            if response == "Y":
                logger.warning("You have chosen to continue, overwriting previous input files.")
            else:
                logger.warning("Continuing with previous input, will not rerun")
                skip = True
    else:
        # Make sure all subdirs exist
        input_path.touch(exist_ok=True)
    with open(input_path, 'w') as f:
        f.write(tmp)
    return input_path, skip, Om_list, w0_list



def reference_setup(options, config, logger):
    if options is None:
        logger.error("No reference cosmology specified, quitting")
        return None
    name = list(options.keys())[0]
    options = options[name]
    Om = options["OMEGA_MATTER"]
    w0 = options["W0_LAMBDA"]
    input_path = config.get("pippin_inputs") / f"BV_FR_{name}.yml"
    master = config.get("input_files") / "reference_master.yml"
    Om_min = Om - 1
    Om_max = Om + 1
    w0_min = w0 - 2
    w0_max = w0 + 2
    # Edit the master yml file
    with open(master, 'r') as f:
        tmp = f.read()
    prepend = f"# Generated via:\n# [ reference.{name} ]\n# OMEGA_MATTER = {Om}\n# W0_LAMBDA = {w0}"
    tmp = tmp.format(**{"PREPEND": prepend, "W0_LAMBDA": w0, "OMEGA_MATTER": Om, "OMEGA_LAMBDA": 1 - Om, "Om_min": Om_min, "Om_max": Om_max, "w0_min": w0_min, "w0_max": w0_max})
    # Check if a file already exists, and whether you should overwrite that file
    skip = False
    if input_path.exists():
        logger.info(f"Found another input file under {name}, comparing them now")
        with open(input_path, 'r') as f:
            tmp2 = f.read()
        diffs = list(difflib.unified_diff(tmp.split('\n'), tmp2.split('\n'), fromfile="new_input", tofile="old_input", lineterm=''))
        if len(diffs) == 0:
            logger.info("No differences found, will not rerun")
            skip = True
        else:
            logger.warning(f"Found the following differences:")
            for line in diffs:
                logger.warning(line)
            logger.warning("If you continue, you will overwrite previously defined inputs. Would you like to continue (y/n)?")
            response = input().upper()
            if response == "Y":
                logger.warning("You have chosen to continue, overwriting previous input file.")
            else:
                logger.warning("Continuing with previous input, will not rerun")
                skip = True
    else:
        # Make sure all subdirs exist
        input_path.touch(exist_ok=True)
    with open(input_path, 'w') as f:
        f.write(tmp)
    return input_path, skip, Om, w0

def analyse(pippin_outputs, options, config, logger):
    if len(options.keys()) == 0:
        logger.warning("No analysis options chosen, quitting")
        return None
    #num_cpu = multiprocessing.cpu_count()
    num_cpu = 8
    logger.info(f"Num CPU: {num_cpu}")
    wfits = {
        "FR": {},
        "V": {}
    }
    for (name, (path, Om_l, w0_l)) in pippin_outputs.items():
        if "_FR_" in name:
            logger.info("Loading in full runthrough wfit output")
            d = wfits["FR"]
        else:
            logger.info("Loading in validation wfit output")
            d = wfits["V"]
        wfit_dirs = get_wfit_dirs(path)
        wfit_file_dict = {wfit_dir: get_wfit_files(wfit_dir) for wfit_dir in wfit_dirs}
        for (i, wfit_dir) in enumerate(wfit_file_dict.keys()):
            if "_FR_" in name:
                Om = Om_l
                w0 = w0_l
            else:
                Om = Om_l[i]
                w0 = w0_l[i]
            wfit_files = wfit_file_dict[wfit_dir]
            d[str(wfit_dir)] = {}
            with Pool(num_cpu) as pool:
                for (f, data, best) in tqdm(pool.imap(read_wfit, wfit_files), total=len(wfit_files)):
                    d[str(wfit_dir)][str(f)] = {}
                    d[str(wfit_dir)][str(f)]["data"] = data
                    d[str(wfit_dir)][str(f)]["best"] = best
                    d[str(wfit_dir)][str(f)]["Om"] = Om
                    d[str(wfit_dir)][str(f)]["w0"] = w0

    plot_nominal_options = options.get("plot_nominal", None)
    if plot_nominal_options is not None:
        plot_nominal(wfits, plot_nominal_options, config["plot_output"])
    plot_GPE_options = options.get("plot_GPE", None)
    if plot_GPE_options is not None:
        plot_GPE(wfits, plot_GPE_options, config["plot_output"])
    plot_KDE_options = options.get("plot_KDE", None)
    if plot_KDE_options is not None:
        plot_KDE(wfits, plot_KDE_options, config["plot_output"])

def run():
    args = get_args()
    toml_path = Path(args.toml).resolve()
    with open(toml_path, 'rb') as f:
        toml = tml.load(f)
    # Setup paths
    config["base"] = Path(config["base"])
    config["files"] = config["base"] / "files"
    config["base_files"] = config["files"] / "base_files"
    config["input_files"] = config["files"] / "input_files"
    config["inputs"] = config["base"] / "inputs"
    config["pippin_inputs"] = config["inputs"] / "Pippin"
    config["outputs"] = config["base"] / "outputs"
    config["toml_input"] = toml_path.parents
    config["output"] = config["outputs"] / toml_path.stem
    config["plot_output"] = config["output"] / "plots"
    config["pippin_output"] = Path(os.environ["PIPPIN_OUTPUT"])
    config["pippin_outputs"] = {}
    output = config["output"]
    if not output.exists():
        output.mkdir(parents=True, exist_ok=True)
    if not config["plot_output"].exists():
        config["plot_output"].mkdir(parents=True, exist_ok=True)
    logfile = output / config.get("logfile")
    if not logfile.exists():
        logfile.touch(exist_ok=True)
    setup_logging(logfile, args.verbose)
    logger = logging.getLogger("bias_validation")
    logger.info(f"Config complete, output folder set to {output}, logging to {logfile})")
    logger.info("Preparing reference cosmology pippin input")

    # Start of reference cosmology
    reference_options = toml.get("reference", None)
    if reference_options is not None:
        reference_path, skip_reference, reference_Om, reference_w0 = reference_setup(reference_options, config, logger)
        if reference_path is None:
            logger.warning("No full runthrough paths defined")
        else:
            if skip_reference and not args.no_skip:
                logging.info(f"Skipping {reference_path}")
            else:
                logging.info(f"Reference setup complete, running {reference_path}")
                cmd = ["pippin.sh", "-v", str(reference_path)]
                logger.debug(f"Running {cmd}")
                rtn = subprocess.run(cmd).returncode
                
                logger.debug(f"Completed with exit code {rtn}")
                if rtn != 0:
                    logger.warning(f"Non-zero return code {rtn}. The reference job likely failed, will attempt to continue anyway")
            job_name = reference_path.stem.upper()
            output_dir = config["pippin_output"] / job_name
            if output_dir.exists():
                config["pippin_outputs"][job_name] = (output_dir, reference_Om, reference_w0)
                logger.debug(f"Added {output_dir} to analysis list")
            logger.info("Completed reference cosmology")

    # Validation jobs
    validation_options = toml.get("validation", None)
    if validation_options is not None:
        validation_path, skip_validation, validation_Om, validation_w0 = validation_setup(validation_options, config, logger) 
        if validation_path is None:
            logger.warning("No validation paths defined")
        else:
            if skip_validation and not args.no_skip:
                logging.info(f"Skipping {validation_path}")
            else:
                logging.info(f"validation setup complete, running {validation_path}")
                cmd = ["pippin.sh", "-v", str(validation_path)]
                logger.debug(f"Running {cmd}")
                rtn = subprocess.run(cmd).returncode
                
                logger.debug(f"Completed with exit code {rtn}")
                if rtn != 0:
                    logger.warning(f"Non-zero return code {rtn}. The validation job likely failed, will attempt to continue anyway")
            job_name = validation_path.stem.upper()
            output_dir = config["pippin_output"] / job_name
            if output_dir.exists():
                config["pippin_outputs"][job_name] = (output_dir, validation_Om, validation_w0)
                logger.debug(f"Added {output_dir} to analysis list")
            logger.info("Completed validation cosmology")

    # Analysis
    if toml.get("analyse", None) is not None:
        logger.debug(config["pippin_outputs"])
        options = toml["analyse"]
        analyse(config["pippin_outputs"], options, config, logger)

if __name__ == "__main__":
    run()
