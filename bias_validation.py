#!/usr/bin/env python
import argparse
import logging
import coloredlogs
import tomli as tml
from pathlib import Path
import difflib
import shutil
import subprocess
import sys

# Relative imports
from files.scripts.confidence_intervals import calculate_confidence_interval, analyse_confidence_interval
from files.scripts.contour import get_chain_paths, load_chains, plot_chains

config = {
        'logfile': 'log.txt',
        'base': '/project2/rkessler/SURVEYS/DES/USERS/parmstrong/inputs/bias_validation/',
        'clean': True
        }


def ensure_list(l):
    if isinstance(l, list):
        return l
    return [l]

# Logging stuff just stolen from pippin

def setup_logging(logging_filename, verbose):

    level = logging.DEBUG if verbose else logging.INFO

    fmt = "[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s" if verbose else "%(message)s"
    handlers = [logging.StreamHandler()]
    handlers[0].setLevel(level)
    handlers.append(logging.FileHandler(logging_filename, mode="w"))
    handlers[-1].setLevel(logging.DEBUG)
    logging.basicConfig(level=level, format=fmt, handlers=handlers)

    coloredlogs.install(
        level=level,
        fmt=fmt,
        reconfigure=True,
        level_styles=coloredlogs.parse_encoded_styles("debug=8;info=green;warning=yellow;error=red,bold;critical=red,inverse"),
    )
    logging.getLogger("matplotlib").setLevel(logging.ERROR)


def full_run_setup(full_runs, config, logger):
    if len(full_runs) == 0:
        logger.error("No full runthroughs specified, quitting")
        return None
    full_run_paths = []
    for full_run in full_runs:
        logger.info(f"Setting up {full_run}")
        # Grab parameters from input file
        Om = full_runs[full_run][0]['OMEGA_MATTER']
        w0 = full_runs[full_run][0]['W0_LAMBDA']
        name = f"BV_FR_{full_run}.yml"
        full_run_master = config.get("input_files") / "full_run.yml"
        # Edit the master yml file
        with open(full_run_master, 'r') as f:
            tmp = f.read()
        tmp =tmp.replace("$W0_LAMBDA", str(w0))
        tmp = tmp.replace("$OMEGA_MATTER", str(Om))
        # Write out to pippin_inputs
        prepend = f"#Generated via:\n#[[ full_run.{full_run} ]]\n#OMEGA_MATTER = {Om}\n#W0_LAMBDA = {w0}\n\n"
        tmp = prepend + tmp
        full_run_final = config.get("pippin_inputs") / name
        skip = False # Whether or not to skip running pippin
        if full_run_final.exists():
            logger.info(f"Found another input file called {name}, comparing them now")
            with open(full_run_final, 'r') as f:
                tmp2 = f.read()
            diffs = list(difflib.unified_diff(tmp.split('\n'), tmp2.split('\n'), fromfile="new_input", tofile="old_input", lineterm=''))
            if len(diffs) == 0:
                logger.info("No differences found, will not rerun")
                skip = True
            else:
                logger.warning(f"Found the following differences:")
                for line in diffs:
                    logger.warning(line)
                logger.warning("If you continue, you will overwrite previously defined inputs. Would you like to continue? (y/n)")
                response = input()
                if response == 'y':
                    logger.warning("You have chosen to continue, overwriting previous input file.") 
                else:
                    logger.warning("Quitting")
                    return None 
        if not full_run_final.exists():
            full_run_final.touch(exist_ok=True)
        with open(full_run_final, 'w') as f:
            f.write(tmp)
        full_run_paths.append((full_run_final, skip, Om, w0))
    return full_run_paths

def validation_setup(validations, config, logger):
    if len(validations) == 0:
        logger.error("No validations specified, quitting")
        return None
    validation_paths = []
    for validation in validations:
        logger.info(f"Setting up {validation}")
        # Grab parameters from input file
        Om = validations[validation][0]["OMEGA_MATTER"]
        w0 = validations[validation][0]["W0_LAMBDA"]
        ompri = str(validations[validation][0].get("ompri", True))
        if ompri in ["True", "true", "t", "T", "TRUE"]: # Probably overkill
            ompri = "True"
        else:
            try:
                ompri = float(ompri)
                if (ompri < 0) or (ompri > 1):
                    logger.error(f"ompri values: {ompri} must be between 0 and 1")
                    return None
                ompri = str(ompri)
            except ValueError as e:
                logger.error(f"Unknown ompri option: {ompri}. Possible values are True (to set ompri to the simulated omega matter) or a value between 0 and 1")
                return None
        dompri = str(validations[validation][0].get("dompri", 0.05))
        name = f"BV_V_{validation}.yml"
        validation_master = config.get("input_files") / "validation.yml"
        SALT2mu_master = config.get("base_files") / "SALT2mu_master.input"
        DESSIM = ""
        BCOR = ""
        for m in ensure_list(Om):
            m_str = str(m).replace('-','n').replace('.','_')
            min_w = min(ensure_list(w0)) - 1
            max_w = max(ensure_list(w0)) + 1
            WFIT = f"-ompri {str(m) if ompri == 'True' else ompri} -dompri {dompri} -wmin {min_w} -wmax {max_w} -wsteps 201 -hsteps 121"
            with open(SALT2mu_master, 'r') as f:
                tmp = f.read()
            tmp = tmp.replace('$WFIT', WFIT)
            min_w = str(min_w).replace('-', 'n').replace('.','_')
            max_w = str(max_w).replace('-', 'n').replace('.','_')
            with open(config.get("toml_input") / f"SALT2mu_m_{m_str}_w_{min_w}_{max_w}.input", 'w') as f:
                f.write(tmp)
            for w in ensure_list(w0):
                w_str = str(w).replace('-','n').replace('.','_')
                s = f"  DESSIM_w_{w_str}_om_{m_str}:\n    <<: *sim_config\n    GLOBAL:\n      <<: *sim_global\n      OMEGA_MATTER: {str(float(m))}\n      W0_LAMBDA: {str(float(w))}\n      OMEGA_LAMBDA: {str(1 - float(m))}\n\n"
                DESSIM += s
                s = f"  BCOR_w_{w_str}_om_{m_str}:\n      <<: *bcor_config\n      BASE: SALT2mu_m_{m_str}_w_{min_w}_{max_w}.input\n      DATA: [DESFITSYS_DESSIM_w_{w_str}_om_{m_str}]\n\n"
                BCOR += s
        with open(validation_master, 'r') as f:
            tmp = f.read()
        tmp = tmp.replace('$DESSIM', DESSIM)
        tmp = tmp.replace('$BCOR', BCOR)
        tmp = tmp.replace('$DATA_DIR', str(config.get("toml_input")))
        # Write out to pippin_inputs
        prepend = f"#Generated via:\n#[[ validation.{validation} ]]\n#ompri = {ompri}\n#dompri = {dompri}\n#OMEGA_MATTER = {Om}\n#W0_LAMBDA = {w0}\n\n"
        tmp = prepend + tmp
        validation_final = config.get("pippin_inputs") / name
        skip = False # Whether or not to skip running pippin
        if validation_final.exists():
            logger.info(f"Found another input file called {name}, comparing them now")
            with open(validation_final, 'r') as f:
                tmp2 = f.read()
            diffs = list(difflib.unified_diff(tmp.split('\n'), tmp2.split('\n'), fromfile="new_input", tofile="old_input", lineterm=''))
            if len(diffs) == 0:
                logger.info("No differences found, will not rerun")
                skip = True
            else:
                logger.warning(f"Found the following differences:")
                for line in diffs:
                    logger.warning(line)
                logger.warning("If you continue, you will overwrite previously defined inputs. Would you like to continue? (y/n)")
                response = input()
                if response == 'y':
                    logger.warning("You have chosen to continue, overwriting previous input file.")
                else:
                    logger.warning("Quitting.")
                    return None
        if not validation_final.exists():
            validation_final.touch(exist_ok=True)
        with open(validation_final, 'w') as f:
            f.write(tmp)
        validation_paths.append((validation_final, skip))
    return validation_paths
    
def get_args():
    # Set up command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("toml", help="the name of the toml config file to run.", type=str)
    parser.add_argument("-v", "--verbose", help="Increase output verbosity", action="store_true")
    parser.add_argument("-n", "--no_skip", help="Whether to rerun pippin jobs regardless of whether a change has been made", action="store_true")

    args = parser.parse_args()
    return args

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
    config["toml_input"] = toml_path.parent
    config["output"] = config["outputs"] / toml_path.stem # Add config file name to output
    output = config["output"]
    if not output.exists():
        output.mkdir(parents=True, exist_ok=True)
    logfile = output / config.get("logfile")
    if not logfile.exists():
        logfile.touch(exist_ok=True)
    setup_logging(logfile, args.verbose)
    logger = logging.getLogger("bias_validation") 
    logger.info(f"Config complete, output folder set to {output}, logging to {logfile})")
    logger.info("Preparing for full run stage")
    # Start of run
    full_run_paths = full_run_setup(toml.get("full_run", []), config, logger)
    if full_run_paths is None:
        logger.error("Error occured during full runthrough setup, quitting.")
        return None
    logger.info(f"Full runthrough setup complete, found {len(full_run_paths)} runthroughs")
    logger.debug(f"Full runthrough paths: {full_run_paths}")
    logger.info("Running the full runthroughs now")
    for i, (path, skip, _, _) in enumerate(full_run_paths):
        if skip and not args.no_skip:
            logger.info(f"Skipping {i + 1}")
            continue
        logger.info(f"({i + 1}/{len(full_run_paths)})")
        cmd = ["pippin.sh", str(path)]
        logger.debug(f"Running {cmd}")
        cmd = subprocess.run(cmd)
        logger.debug(f"Completed with exit code {cmd.returncode}")
    logger.info("Completed full runthroughs")
    # Start of validation
    logger.info("Setting up validation jobs")
    validation_paths = validation_setup(toml.get("validation", []), config, logger)
    if validation_paths is None:
        logger.error("Error occured during validation jobs setup, quitting")
        return None
    logger.info(f"Validation setup complete, found {len(validation_paths)} jobs")
    logger.debug(f"Validation jobs: {validation_paths}")
    logger.info("Running validation jobs now")
    for i, (path, skip) in enumerate(validation_paths):
        if skip and not args.no_skip:
            logger.info(f"Skipping {i + 1}")
            continue
        logger.info(f"({i + 1}/{len(validation_paths)})")
        cmd = ["pippin.sh", str(path)]
        logger.debug(f"Running {cmd}")
        cmd = subprocess.run(cmd)
        logger.debug(f"Compeleted with exit code {cmd.returncode}")
    logger.info("Completed validation jobs")
    # Calculating confidence intervals
    logger.info("Calculating confidence intervals")
    ci_config = toml["confidence_interval"]
    do_plot = ci_config.get("plot", False) 
    plot_path = output if do_plot else None
    ci = {}
    for (path, skip) in validation_paths:
        name = path.stem.upper()    
        ci[name] = {}
        ci_stats = calculate_confidence_interval(name, logger)
        for w0 in ensure_list(ci_config["w0"]):
            for Om in ensure_list(ci_config["Om"]):
                wl, wr, oml, omr = analyse_confidence_interval(*ci_stats, name, w0, Om, logger, plot_path)
                ci[name]["w0l"] = wl
                ci[name]["w0r"] = wr
                ci[name]["oml"] = oml
                ci[name]["omr"] = omr
    logger.debug(f"Confidence intervals: {ci}")
    logger.info("Finished calculating confidence intervals")
    # Creating contours 
    logger.info("Creating full runthrough contours")
    for (path, skip, Om, w0) in full_run_paths:
        logger.info(f"Creating contour for {path.stem.upper()}")
        chain_paths = get_chain_paths(path.stem.upper(), logger)
        chains = load_chains(chain_paths, masks=[['(SN)'], ['(SN + PRIOR)']], logger=logger) 
        plot_chains(output, chains, w0=w0, Om=Om, confidence_interval=ci, logger=logger)
    logger.info(f"All done! See {output} for output")



if __name__ == "__main__":
    run()
