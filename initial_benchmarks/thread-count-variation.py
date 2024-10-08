#!/usr/bin/env python3

# 
# Script to visualize results produced by process-results.py.
# Created to produce graphs for thread-count variation experiment;
# Also denoted as experiment 2.
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

def create_plot(thread_count: List[int], times: Dict[str, float], options, log: bool=False) -> None:
    """ Creates the default plot for experiment 2: thread count variation. """
    datapath: str = options["datapath"]
    
    x = np.arange(len(thread_count))

    fig, ax = plt.subplots()
    # fig.tight_layout()
    markers: str = "x+*2"
    count: int = 0

    for implementation, time in times.items():
        if len(thread_count) != len(time):
            print(f"ERROR: {implementation} has a mismatch in thread counts ({len(thread_count)}) and time ({len(time)})", file=sys.stderr)
        plt.plot(thread_count, time, marker=markers[count], label=implementation)
        count += 1
    
    if log:
        plt.yscale("log")

    ax.set_ylabel("Calculation time (s)")
    ax.set_xlabel("Thread count")
    # ax.set_title(f"Thread count variation ({options['compute-cluster']})")
    ax.set_xticks(thread_count)
    ax.legend(loc="upper right", ncol=1)
    # ax.set_ylim(0, 3)

    if log:
        plt.savefig(f"{datapath}/ex{options['ex-num']}_{options['compute-cluster']}_tcv_log_{options['file-date']}.png")
    else:
        plt.savefig(f"{datapath}/ex{options['ex-num']}_{options['compute-cluster']}_tcv_{options['file-date']}.png")
    plt.close(fig)
    # plt.show()

def main():
    # Header variable to use as x-axis in result graphs
    # Can be chosen from 'headers' present in common.py
    variable: str = "Thread Count"

    # Argument parsing
    parser: ArgumentParser = ArgumentParser()
    common.add_parameters(parser)
    parser.add_argument('--log', action='store_true')
    args = parser.parse_args()
    options = dict(vars(args))

    prefix = f"EX{options['ex-num']}-{options['compute-cluster']}-RESULTS-"
    
    # Check if results file exists
    file: str = options["datapath"] + "/" + prefix + options["file-date"] + options["file-extension"]
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
    data = common.average_calc_times(data)
    thread_count, times = common.serialize(data)

    # Create Plots
    log: bool = True if args.log else False
    create_plot(thread_count, times, options, log)
    data = common.process_reverse(df, variable)
    common.create_confidence_interval(data, options, 8, options["ex-num"])

if __name__ == "__main__":
    sys.exit(main())