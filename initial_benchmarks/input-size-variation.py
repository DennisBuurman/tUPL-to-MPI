#!/usr/bin/env python3

# 
# Script to visualize results produced by process-results.py.
# Created to produce graphs for input size variation experiment;
# Also denoted as experiment 1.
# 
# Author: Dennis Buurman, Leiden University

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import numpy as np
from argparse import ArgumentParser

import common

from typing import Dict, List

def create_plot(sizes: List[int], times: Dict[str, float], options) -> None:
    """ Creates the default plot for experiment 1: input size variation. """
    datapath: str = options["datapath"]
    ex_num: int = options["ex-num"]

    x = np.arange(len(sizes))
    width: float = 0.2 # bar width
    multiplier: int = 0

    fig, ax = plt.subplots()
    fig.tight_layout()

    for implementation, time in times.items():
        offset: float = width * multiplier
        ax.bar(x + offset, time, width, label=implementation)
        # ax.bar_label(rects, padding=3)
        multiplier += 1
    
    ax.set_ylabel("Calculation time (s)")
    ax.set_xlabel("Input size (2^x)")
    ax.set_title(f"Input size variation: {options['compute-cluster']}")
    ax.set_xticks(x + width)
    ax.set_xticklabels(sizes)
    ax.legend(loc="upper left", ncol=1)
    # ax.set_ylim(0, 3)

    plt.savefig(f"{datapath}/ex{ex_num}_{options['compute-cluster']}_{options['file-date']}.png")
    plt.close(fig)
    # plt.show()

def create_boxplots(data: Dict, options) -> None:
    """ Creates boxplots for each implementation to visualize performance variability. """
    cluster: str = options["compute-cluster"]
    date: str = options["file-date"]
    datapath: str = options["datapath"]
    ex_num: int = options["ex-num"]

    for implementation in data:
        fig, ax = plt.subplots()
        fig.tight_layout()
        box_plot_data = []
        labels = []
        for size in data[implementation]:
            times = data[implementation][size]
            box_plot_data.append(times)
            labels.append(size)
        ax.set_title(f"{common.names[implementation]} {cluster} runtime distribution 2^{size}")
        ax.set_ylabel("Calculation time (s)")
        ax.set_xlabel("Input size (2^x)")
        ax.boxplot(box_plot_data, patch_artist=True, labels=labels)
        plt.savefig(f"{datapath}/ex{ex_num}_bp_{implementation}_{cluster}_{date}")
        plt.close(fig)

def main():
    variable: str = "Input Size"

    # Argument parsing
    parser: ArgumentParser = ArgumentParser()
    common.add_parameters(parser)
    args = parser.parse_args()
    options: Dict[str, any] = dict(vars(args))
    
    prefix = f"EX{options['ex-num']}-{options['compute-cluster']}-RESULTS-"

    file: str = options["datapath"] + "/" + prefix + options["file-date"] + options["file-extension"]
    if not Path(file).is_file():
        print(f"Error: {file} is not a file.", file=sys.stderr)
        return 1
    
    # Preprocessing
    df: DataFrame = pd.read_csv(file, delim_whitespace=True, names=common.headers)
    data: Dict[str, Dict[str, List[float]]] = common.process(df, variable)
    data = common.average_calc_times(data)
    sizes, times = common.serialize(data)

    # Create plots
    create_plot(sizes, times, options)
    data = common.process_reverse(df, variable)
    create_boxplots(data, options)
    common.create_confidence_interval(data, options, 28, options['ex-num'])

    return 0

if __name__ == "__main__":
    sys.exit(main())