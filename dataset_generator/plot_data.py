#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser

from typing import Dict, List, Tuple

datafile = "data.txt"
meansfile = "initial_means{}.txt"


def plot2d(d: np.array, m: np.array, meta: Dict[str, int]) -> None:
    """ Create 2d plot of dataset """
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    # Add data points
    x = d[:, 0]
    y = d[:, 1]
    ax1.scatter(x, y, alpha=0.7, c='b')

    # Color initial means
    x = m[:, 0]
    y = m[:, 1]
    ax1.scatter(x, y, alpha=1.0, c='r', label="Initial means")
    plt.title("k-means 2d scatter plot")
    plt.legend(loc='upper left')
    # plt.show()
    plt.savefig(f"k-means_2d_seed{meta['seed']}_size{meta['size']}_c{meta['clusters']}.png")

def plot3d(d: np.array, m: np.array, meta: Dict[str, int]) -> None:
    print("Not implemented yet")
    pass

def plot4d(d: np.array, m: np.array, meta: Dict[str, int]) -> None:
    print("Not implemented yet")
    pass

def plot_dataset(directory: str, means: int, meta: Dict[str, int]) -> None:
    """ Plot the data from datafilename and highlight the initial means from meansfilename. """
    with open(directory+datafile, "r") as f:
        d = np.array([[float(y) for y in x.split()[1:]] for x in f.readlines()])
    with open(directory+meansfile.format(means), "r") as f:
        m = np.array([[float(y) for y in x.split()[1:]] for x in f.readlines()])
    
    dim = d.shape[1]
    plotting_functions: Dict[int, function] = {
        2: plot2d,
        3: plot3d,
        4: plot4d
    }
    if (dim not in plotting_functions):
        print(f"ERROR: unable to plot dataset of dimension {dim}", file=sys.stderr)
    plotting_functions[dim](d, m, meta)

def read_meta(dir: str) -> Dict[str, int]:
    """ Read the metadata file from the given directory """
    meta = {}
    with open(dir+"metadata.md", "r") as f:
        for x in f.readlines()[1:]:
            split = x.split(":")
            meta[split[0]] = int(split[1])
    return meta

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-d", "--directory", dest="directory", type=str, required=True,
                        help="Path to dataset directory")
    parser.add_argument("-m", "--means_set", dest="means_set", type=int, required=True,
                        help="Path to dataset directory")
    args = parser.parse_args()
    options: Dict[str, any] = dict(vars(args))

    meta = read_meta(options["directory"])
    plot_dataset(options["directory"], options["means_set"], meta)
