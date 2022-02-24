#!/usr/bin/env python

import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from pathlib import Path
import argparse

def calculate_confidence_interval(name, logger=None):
    # name: Name of pippin folder in $PIPPIN_OUTPUT
    # logger: Optional logging
    # plot: If None, don't plot, otherwise, should equal the output directory the plots will go in.
    #       The plots will be named "{name}_w0.png" and "{name}_Om.png"


    ROOT_DIR = Path(f"/scratch/midway2/rkessler/PIPPIN_OUTPUT/{name}/6_BIASCOR")

    cosmology_dirs = [p for p in ROOT_DIR.iterdir() if p.is_dir()]

    realisations_with_missing_cospar = []

    w_true = []
    om_true = []
        
    w_obs_ci_min = []
    w_obs_ci_max = []
    om_obs_ci_min = []
    om_obs_ci_max = []
    w_dtype = [('w_true', float), ('w_obs_ci_min', float), ('w_obs_ci_max', float)]
    om_dtype = [('om_true', float), ('om_obs_ci_min', float), ('om_obs_ci_max', float)]
    # WIP: structured array to keep w_obs and w_true values mapped to each other

    for cosmology_dir in cosmology_dirs:
        cosmology_dir_path = cosmology_dir / "output"
        subdirs_to_process = [p for p in cosmology_dir_path.iterdir() if p.is_dir and "BV_V" in p.name]

        wvals = []
        werr = []
        omvals = []
        omerr = []

        for subdir in subdirs_to_process:
            yaml_path = subdir / "wfit_FITOPT000_MUOPT000.YAML"

            if not yaml_path.is_file():
                realisations_with_missing_cospar.append(subdir.stem)
                continue

            # Read w and om values from file
            with open(yaml_path) as f:
                data = f.readlines()
                w = float(data[0].strip().split(":")[-1])
                w_sig = float(data[1].strip().split(":")[-1])
                om = float(data[2].strip().split(":")[-1])
                om_sig = float(data[3].strip().split(":")[-1])
                wvals.append(w)
                werr.append(w_sig)
                omvals.append(om)
                omerr.append(om_sig)

        if len(wvals) == 0:
            if logger is not None:
                logger.debug("No values found")
            continue

        n = cosmology_dir.name
        n = n[7:] # Remove BCOR_ from name
        w = float(n[:n.index("_om_")].replace('n','-').replace('_','.')) # Recover the value of w from the name
        n = n[n.index("_om_")+4:] # Remove w info and _om_ from name
        om = float(n.replace("n","-").replace("_","."))
        w_true.append(w)
        om_true.append(om)

        # Calculate w statistics
        w_weight = 1/(np.array(werr)**2)
        w_avg = np.average(wvals,weights=w_weight)
        w_err = np.average(werr,weights=w_weight)
        w_var = np.sqrt(1/np.sum(w_weight))
        w_std = np.std(wvals)

        w_obs_ci_min.append(w_avg - w_err)
        w_obs_ci_max.append(w_avg + w_err)

        # Calculate Om statistics
        om_weight = 1 / (np.array(omerr)**2)
        om_avg = np.average(omvals, weights=om_weight)
        om_err = np.average(omerr, weights=om_weight)
        om_var = np.sqrt(1/np.sum(om_weight))
        om_std = np.std(omvals)

        om_obs_ci_min.append(om_avg - om_err)
        om_obs_ci_max.append(om_avg + om_err)

    if logger is not None and len(realisations_with_missing_cospar) > 0:
        logger.warn(f"Missing realisations in following directories: {realisations_with_missing_cospar}")
        
    # Calculate w 
    w_true = np.array(w_true)
    w_obs_ci_max = np.array(w_obs_ci_max)
    w_obs_ci_min = np.array(w_obs_ci_min)
    mask = np.argsort(w_true)
    w_true_sorted = w_true[mask]
    w_obs_ci_max_sorted = w_obs_ci_max[mask]
    w_obs_ci_min_sorted = w_obs_ci_min[mask]
    
    # Calculate Om
    om_true = np.array(om_true)
    om_obs_ci_max = np.array(om_obs_ci_max)
    om_obs_ci_min = np.array(om_obs_ci_min)
    mask = np.argsort(om_true)
    om_true_sorted = om_true[mask]
    om_obs_ci_max_sorted = om_obs_ci_max[mask]
    om_obs_ci_min_sorted = om_obs_ci_min[mask]

    return w_obs_ci_min_sorted, w_obs_ci_max_sorted, w_true_sorted, om_obs_ci_min_sorted, om_obs_ci_max_sorted, om_true_sorted

def analyse_confidence_interval(w_ci_min, w_ci_max, w_true, om_ci_min, om_ci_max, om_true, name, w0=-1, Om=0.3, logger=None, plot=None):
    w_interp_right = np.interp(w0, w_ci_max, w_true)
    w_interp_left = np.interp(w0, w_ci_min, w_true)

    om_interp_right = np.interp(Om, om_ci_max, om_true)
    om_interp_left = np.interp(Om, om_ci_min, om_true)

    w0_name = str(w0).replace('-','n').replace('.','_')
    Om_name = str(Om).replace('-','n').replace('.','_')

    if plot is not None:
        # w plot
        if len(set(w_true)) > len(set(om_true)):
            plt.hlines(w_true, w_ci_min, w_ci_max)
            plt.axvline(w0, color='r', linestyle='dotted')
            plt.axhline(w_interp_left, color='r', linestyle='dotted')
            plt.axhline(w_interp_right, color='r', linestyle='dotted')
            plt.plot(w_ci_min, w_true, color='red')
            plt.plot(w_ci_max, w_true, color='red')
            plt.title(f"{name}")
            plt.xlabel("Observed w0")
            plt.ylabel("Input w0")
            plt.savefig(plot / f"{name}_w0_{w0_name}_ladder.png")
            plt.clf()

            plt.vlines(w_true, w_ci_min - w_true, w_ci_max - w_true)
            plt.plot(w_true, w0 - w_true, color='red', linestyle='dotted')
            plt.axvline(w_interp_left, color='red', linestyle='dotted')
            plt.axvline(w_interp_right, color='red', linestyle='dotted')
            plt.plot(w_true, w_ci_min - w_true, color='red')
            plt.plot(w_true, w_ci_max - w_true, color='red')
            plt.title(f"{name}")
            plt.xlabel("Input w0")
            plt.ylabel("Observed w0 - input w0")
            plt.savefig(plot / f"{name}_w0_{w0_name}_residual.png")
            plt.clf()

        # Om plot
        if len(set(om_true)) > len(set(w_true)):
            plt.hlines(om_true, om_ci_min, om_ci_max)
            plt.axvline(Om, color='r', linestyle='dotted')
            plt.axhline(om_interp_left, color='r', linestyle='dotted')
            plt.axhline(om_interp_right, color='r', linestyle='dotted')
            plt.plot(om_ci_min, om_true, color='red')
            plt.plot(om_ci_max, om_true, color='red')
            plt.title(f"{name}")
            plt.xlabel("Observed Om")
            plt.ylabel("Input Om")
            plt.savefig(plot / f"{name}_Om_{Om_name}_ladder.png")
            plt.clf()

            plt.vlines(om_true, om_ci_min - om_true, om_ci_max - om_true)
            plt.plot(om_true, Om-om_true, color='red', linestyle='dotted')
            plt.axvline(om_interp_left, color='red', linestyle='dotted')
            plt.axvline(om_interp_right, color='red', linestyle='dotted')
            plt.plot(om_true, om_ci_min - om_true, color='red')
            plt.plot(om_true, om_ci_max - om_true, color='red')
            plt.title(f"{name}")
            plt.xlabel("Input Om")
            plt.ylabel("Observed Om - input Om")
            plt.savefig(plot / f"{name}_Om_{Om_name}_residual.png")
            plt.clf()

    return w_interp_left, w_interp_right, om_interp_left, om_interp_right

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="The name of the pippin output folder", type=str)
    parser.add_argument("-p", "--plot", help="Directory to store plots", type=str, default=None)
    parser.add_argument("-w", "--w0", help="w0 value of interest", type=float, default=-1.0)
    parser.add_argument("-m", "--Om", help="Om value of interest", type=float, defualt=0.3)

    args = parser.parse_args()
    return args

def main():
    args = get_args()
    name = args.name
    plot = args.plot
    Om = args.Om
    w0 = args.w0
    stats = confidence_interval(name)
    wl, wr, oml, omr = analyse_confidence_interval(*stats, name, w0, Om, None, plot)

if __name__ == "__main__":
    main()

