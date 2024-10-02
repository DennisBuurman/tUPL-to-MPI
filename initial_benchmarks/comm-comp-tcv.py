#!/usr/bin/env python3

# 
# Script to visualize results produced by process-results.py.
# Created to produce graphs for communication vs computation thread-count variation experiment;
# Also denoted as experiment 9.
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

from typing import Dict, List, Tuple

variables = [
    "Worst Comm. Time", 
    "Worst Reassign Time"
]

def average_results(data: Dict[str, List[List[float]]]) -> Dict[str, List[List[float]]]:
    """ Average all communication and computation times of all implementations.
        The implementation will be named 'average'. """
    d: Dict[str, List[List[float]]] = {"average": [[], []]} # [[comp], [comm]]

    for implementation in data:
        comp_times = data[implementation][0]
        comm_times = data[implementation][1]
        # first implementation is copied directly
        if len(d["average"][0]) == 0:
            d["average"][0] = comp_times
            d["average"][1] = comm_times
        # rest is summed
        else:
            for i in range(len(comp_times)):
                d["average"][0][i] += comp_times[i]
                d["average"][1][i] += comm_times[i]
    
    count = len(data.keys())
    for i in range(len(d["average"][0])):
        d["average"][0][i] = float(d["average"][0][i]/count)
        d["average"][1][i] = float(d["average"][1][i]/count)
    
    return d
    

def preprocess(df: DataFrame, variable: str, per_iteration=False, average=False) -> Dict[str, Dict[str, List[List[float]]]]:
    """ Processes data from a result file into a dictionary.
        Preps data for general graph. 
        variable should be in the timing_headers list from common.py."""
    d: Dict[str, Dict[str, List[List[float]]]] = dict()
    df.reset_index()

    for v in variables:
        if v not in list(df.columns):
            print(f"ERROR: variable '{variable}' not in dataframe!", file=sys.stderr)
            return {}
    if variable not in list(df.columns):
        print(f"ERROR: variable '{variable}' not in dataframe!", file=sys.stderr)
        return {}
    
    # Sort data on variable, implementation, worst calc. time
    for _, row in df.iterrows():
        implementation: str = row["Implementation"]
        var: str = row[variable]
        comp_time: float = row["Worst Reassign Time"] if not per_iteration else row["Worst Reassign Time"] / row["Iterations"]
        comm_time: float = row["Worst Comm. Time"] if not per_iteration else row["Worst Comm. Time"] / row["Iterations"]
        if var in d:
            time_per_implementation: Dict[str, List[List[float]]] = d[var]
            if implementation in time_per_implementation:
                time_per_implementation[implementation][0].append(comp_time)
                time_per_implementation[implementation][1].append(comm_time)
            else:
                time_per_implementation[implementation] = [[comp_time], [comm_time]]
        else:
            d[var] = {
                implementation: [
                    [comp_time],
                    [comm_time]
                ]
            }
    
    if average:
        d = average_timings(d)

    return d

def average_times(d: Dict[str, Dict[str, List[List[float]]]]) -> Dict[str, Dict[str, List[float]]]:
    """ Calculates the mean of the times present in a nested dictionary created by preprocess(). 
        List of float values is replaced by its mean. """
    for x in d:
        implementations: Dict[str, List[float]] = d[x]
        for i in implementations:
           comp_times: List[float] = implementations[i][0]
           implementations[i][0] = np.mean(comp_times)
           comm_times: List[float] = implementations[i][1]
           implementations[i][1] = np.mean(comm_times)
    
    return d

def serialize(data: Dict[str, Dict[str, float]]) -> Tuple[List[int], Dict[str, List[float]]]:
    """ Serializes the data for plot input. """
    var_list: List[int] = [] # list of variable values
    times: Dict[str, List[List[float]]] = {}

    for var in data:
        var_list.append(var) # thread counts
        for implementation in data[var]:
            comp_t: float = data[var][implementation][0]
            comm_t: float = data[var][implementation][1]
            n: str = common.names[implementation] # Anne's naming convention
            if n in times:
                times[n][0].append(comp_t)
                times[n][1].append(comm_t)
            else:
                times[n] = [[comp_t], [comm_t]]

    return var_list, times

def create_plot(thread_count: List[int], times: Dict[str, List[float]], options, log: bool=False, average=False) -> None:
    """ Creates the default plot for experiment 2: thread count variation. """
    datapath: str = options["datapath"]
    x = np.arange(len(thread_count))

    fig, ax = plt.subplots()
    # fig.tight_layout()
    markers: str = "x+*2.,vspd"
    count: int = 0

    for implementation, time in times.items():
        if len(thread_count) != len(time[0]):
            print(f"ERROR: {implementation} has a mismatch in thread counts ({len(thread_count)}) and time ({len(time[0])})", file=sys.stderr)
        plt.plot(thread_count, time[0], marker=markers[count], label=f"{implementation}_comp. time")
        plt.plot(thread_count, time[1], marker=markers[count+1], label=f"{implementation}_comm. time")
        count += 2
    
    if log:
        plt.yscale("log")

    ax.set_ylabel("Time (s)")
    ax.set_xlabel("Thread count")
    ax.set_xticks(thread_count)
    ax.legend(loc="upper right", ncol=1)
    # ax.set_ylim(0, 3)

    avg: str = "" if not average else "_avg"

    if log:
        plt.savefig(f"{datapath}/ex{options['ex-num']}_{options['compute-cluster']}_comp_comm_tcv{avg}_log_{options['file-date']}.png")
    else:
        plt.savefig(f"{datapath}/ex{options['ex-num']}_{options['compute-cluster']}_comp_comm_tcv{avg}_{options['file-date']}.png")
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
    parser.add_argument('--average', action='store_true')
    args = parser.parse_args()
    options = dict(vars(args))

    prefix = f"EX{options['ex-num']}-{options['compute-cluster']}-TIMING-RESULTS-"
    
    # Check if results file exists
    file: str = options["datapath"] + "/" + prefix + options["file-date"] + options["file-extension"]
    if not Path(file).is_file():
        print(f"Error: {file} is not a file.", file=sys.stderr)
        return 1
    
    # Preprocessing
    tc_list: List[int] = []
    df: DataFrame = pd.read_csv(file, delim_whitespace=True, names=common.timing_headers)
    for _, row in df.iterrows():
        tc_list.append(row["Nodes"] * row["Tasks Per Node"])
    df.insert(5, variable, tc_list, True)
    data: Dict[str, Dict[str, List[List[float]]]] = preprocess(df, variable)
    data = average_times(data)
    thread_count, times = serialize(data)
    average: bool = True if args.average else False
    if average:
        times = average_results(times)

    # Create Plot
    log: bool = True if args.log else False
    create_plot(thread_count, times, options, log, average=average)

if __name__ == "__main__":
    sys.exit(main())