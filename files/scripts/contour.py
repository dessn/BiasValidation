#!/usr/bin/env python

from matplotlib import pyplot as plt
import matplotlib
matplotlib.use("Agg")
from chainconsumer import ChainConsumer
import numpy as np
from pathlib import Path
import argparse
import logging
from tqdm import tqdm
import io

def get_chain_paths(name, logger=None):
    if logger is not None:
        logger.info(f"Reading in chains for {name}")
    else:
        print(f"Reading in chains for {name}")        
    path = Path(f"/scratch/midway2/rkessler/PIPPIN_OUTPUT/{name}/9_ANALYSE/ALL_OMW")
    chain_paths = [p for p in path.iterdir() if p.is_file() and p.suffixes == ['.csv', '.gz']]
    if logger is not None:
        logger.info(f"Found {len(chain_paths)} chains")
    else:
        print(f"Found {len(chain_paths)} chains")
    return chain_paths

def load_chains(chain_paths, masks=None, num_cores=8, logger=None):
    if masks is None:
        masks = [""]
    else:
        masks = [m[0] for m in masks]
    weights = {m: [] for m in masks}
    likelihoods = {m: [] for m in masks}
    Om = {m: [] for m in masks}
    w0 = {m: [] for m in masks}
    if logger is not None:
        logger.info(f"Loading chains into chainconsumer")
    else:
        print(f"Loading chains into chainconsumer")
    def _load_chain(path):
        return (path, np.transpose(np.loadtxt(path, dtype=float, skiprows=1, delimiter=","))) # Load in. remove header, transpose
    shared_data = [_load_chain(path) for path in tqdm(chain_paths)]
    for (path, data) in shared_data:
        for key in masks:
            if key in path.name:
                weights[key].append(data[0])
                likelihoods[key].append(data[1])
                Om[key].append(data[2])
                w0[key].append(data[3])
    chains = {}
    for key in masks:
        weights[key] = np.concatenate(weights[key])
        likelihoods[key] = np.concatenate(likelihoods[key])
        Om[key] = np.concatenate(Om[key])
        w0[key] = np.concatenate(w0[key])
        c = ChainConsumer()
        c.add_chain([Om[key], w0[key]], parameters = ["$\Omega_{M}$", "$\omega_{0}$"], weights=weights[key], power=2)
        c.configure(statistics="cumulative")
        if logger is not None:
            logger.info(f"Finished combining chains for mask {key}:\n{c.analysis.get_summary()}")
        else:
            print(f"Finished combining chains for mask {key}:\n{c.analysis.get_summary()}")
        chains[key] = c
    return chains

def plot_chains(plot, chains, w0=-1, Om=0.3, confidence_interval=None, logger=None):
    if confidence_interval is None:
        confidence_interval = {}
    for key in chains.keys():
        c = chains[key]
        plt.ioff()
        fig_width_pt = 4 * 240
        inches_pt = 1.0 / 72.27
        fig_width = fig_width_pt * inches_pt
        fig = c.plotter.plot(display=False, truth=[w0, Om], figsize=(fig_width, fig_width / np.sqrt(2)))
        ax = fig.get_axes()
        for ci in confidence_interval.values():
            ax[2].plot([ci["oml"], ci["omr"]], [ci["w0l"], ci["w0r"]], c='black')
        plt.savefig(Path(plot) / f"contour_{key}.png")
        plt.clf()

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="The name of the pippin output folder", type=str)
    parser.add_argument("-p", "--plot", help="Directory to store plots", type=str, default=None)
    parser.add_argument("-w", "--w0", help="w0 value used in simulation", type=float, default=-1.0)
    parser.add_argument("-o", "--Om", help="Om value used in simulation", type=float, default=0.3)
    parser.add_argument("-c", "--confidence_interval", help="Confidence interval. Can be added multiple times, syntax is -c OmegaMatterBottom OmegaMatterTop W0Left W0Right", action='append', nargs='+', type=float, default=[])
    parser.add_argument("-m", "--masks", help="Mask pippin analysis output. Can be added multiple times. Partial match", action='append', nargs='+', type=str, default=None)
    parser.add_argument("-n", "--num_cores", help="How many cores to use. Warning, using to many can get you kicked off midway!", type=int, default=8)

    args = parser.parse_args()
    return args

def main():
    args = get_args()
    name = args.name
    plot = args.plot
    Om = args.Om
    w0 = args.w0
    confidence_interval = {}
    ci = args.confidence_interval
    for i, c in enumerate(ci):
        confidence_interval[i] = {'oml': c[0], 'omr': c[1], 'w0l': c[2], 'w0r': c[3]}
    masks = args.masks
    num_cores = args.num_cores
    chain_paths = get_chain_paths(name)
    chains = load_chains(chain_paths, masks=masks, num_cores=num_cores)
    if plot is not None:
        plot_chains(plot, chains, w0=w0, Om=Om, confidence_interval=confidence_interval)

if __name__ == "__main__":
    main()
