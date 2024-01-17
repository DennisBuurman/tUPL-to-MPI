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

def create_plot(thread_count: List[int], times: Dict[str, float], options) -> None:
    """ Creates the default plot for experiment 2: thread count variation. """
    name: str = f"{options['compute-cluster']}_{options['file-date']}"
    datapath: str = options["datapath"]
    
    x = np.arange(len(thread_count))

    fig, ax = plt.subplots()
    fig.tight_layout()
    markers: str = "x+*2"
    count: int = 0

    for implementation, time in times.items():
        plt.plot(thread_count, time, marker=markers[count], label=implementation)
        count += 1
    
    # plt.yscale("log")

    ax.set_ylabel("Calculation time (s)")
    ax.set_xlabel("Thread count")
    ax.set_title(f"Thread count variation: {name}")
    ax.set_xticks(thread_count)
    ax.legend(loc="upper right", ncol=1)
    # ax.set_ylim(0, 3)

    plt.savefig(f"{datapath}/exp_2_{name}.png")
    plt.close(fig)
    # plt.show()

def main():
    variable: str = "Thread Count"

    # Argument parsing
    parser = ArgumentParser()
    common.add_argument(parser)
    args = parser.parse_args()
    options = dict(vars(args))

    file: str = options["datapath"] + "/" + options["file-preamble"] + options["file-date"] + options["file-extension"]
    if not Path(file).is_file():
        print(f"Error: {file} is not a file.", file=sys.stderr)
        return 1
    
    # Preprocessing
    tc_list: List[int] = []
    df: DataFrame = pd.read_csv(file, delim_whitespace=True, names=common.headers)
    for _, row in df.iterrows():
        tc_list.append(row["Nodes"] * row["Tasks Per Node"])
    df.insert(5, variable, tc_list, True)
    data: Dict = common.process(df, variable)
    thread_count, times = common.serialize(data)

    # Create Plots
    create_plot(thread_count, times, options)
    data = common.process_reverse(df, variable)
    common.create_confidence_interval(data, options, 8)

if __name__ == "__main__":
    sys.exit(main())