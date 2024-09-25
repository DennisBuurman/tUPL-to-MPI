#!/usr/bin/env python3

# 
# Script to visualize results produced by process-timing.py.
# Produces a graph which compares the scaling factor of communication and computation per thread-count.
# Also denoted as experiment 3.
# 
# Author: Dennis Buurman, Leiden University

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import numpy as np
from argparse import ArgumentParser

from pprint import pprint

import common

from typing import Dict, List

# Variables used for plot
# Variables must be in 'timing_headers' from common.py!
variables = [
    "Worst Comm. Time", 
    "Worst Reassign Time"
]

def average_timings(data: Dict[str, Dict[str, Dict[int, List[float]]]]) -> Dict[str, Dict[str, Dict[int, List[float]]]]:
    """ Average all communication and computation times of all implementations.
        The implementation will be named 'average'. """
    d: Dict[str, Dict[str, Dict[int, List[float]]]] = {"average":{}}

    # Sum all comp and comm times from each implementation per thread count
    for implementation in data:
        for var in data[implementation]:
            if var not in d["average"]:
                d["average"][var] = {}
            for tc in data[implementation][var]:
                comp_time: float = data[implementation][var][tc][0]
                comm_time: float = data[implementation][var][tc][1]
                if tc not in d["average"][var]:
                    d["average"][var][tc] = [comp_time, comm_time]
                else:
                    d["average"][var][tc][0] += comp_time
                    d["average"][var][tc][1] += comm_time
    
    # Average all comp and comm times by dividing them with the amount of implementations
    count: int = len(data.keys())
    for var in d["average"]:
        for tc in d["average"][var]:
            d["average"][var][tc][0] = float(d["average"][var][tc][0]/count)
            d["average"][var][tc][1] = float(d["average"][var][tc][1]/count)

    return d

def preprocess(df: DataFrame, average: bool = False) -> Dict[str, Dict[str, List[List[float]]]]:
    """ Preprocessing for scaling plot. Creates a dictionary containing the intra-node and 
        inter-node scaling factors for both worst. computation time and worst. communication 
        time for each implementation. """
    d: Dict[str, Dict[str, List[List[float]]]] = dict() # intermediate typing changes nested list structure to Dict[int, List[float]]
    df.reset_index()

    # Check if variables are present in df
    for v in variables:
        if v not in list(df.columns):
            print(f"ERROR: variable '{v}' not in dataframe!", file=sys.stderr)
            return d
    
    # Preprocess
    for _, row in df.iterrows():
        comm_time: float = row["Worst Comm. Time"]
        comp_time: float = row["Worst Reassign Time"]
        nodes: int = row["Nodes"]
        tasks: int = row["Tasks Per Node"]
        thread_count: int = nodes * tasks
        implementation: str = row["Implementation"]

        # Sort per implementation
        if implementation not in d:
            d[implementation] = dict()
        
        # Determine if intra or inter node (or both)
        if nodes == 1:
            # Intra node
            if "intra" not in d[implementation]:
                d[implementation]["intra"] = {}
            if thread_count not in d[implementation]["intra"]:
                d[implementation]["intra"][thread_count] = [[], []] # 0: comp_time, 1: comm_time
            d[implementation]["intra"][thread_count][0].append(comp_time)
            d[implementation]["intra"][thread_count][1].append(comm_time)  
            if tasks == 1:
                print(f"WARNING: including thread count 1 results in distorted communication time!", file=sys.stderr)
                # Both intra and inter
                if "inter" not in d[implementation]:
                    d[implementation]["inter"] = {} 
                if thread_count not in d[implementation]["inter"]:
                    d[implementation]["inter"][thread_count] = [[], []] # 0: comp_time, 1: comm_time
                d[implementation]["inter"][thread_count][0].append(comp_time)
                d[implementation]["inter"][thread_count][1].append(comm_time) 
        elif nodes > 1:
            # Inter node
            if tasks > 1:
                # Higher tasks per node distorts inter node effect of adding threads
                # Only the case for strictly distributed-memory execution
                print(f"WARNING: Inter node thread counts use >1 '{tasks}' tasks per node!", file=sys.stderr)
            if "inter" not in d[implementation]:
                d[implementation]["inter"] = {} 
            if thread_count not in d[implementation]["inter"]:
                d[implementation]["inter"][thread_count] = [[], []]# 0: comp_time, 1: comm_time, 2: threads
            d[implementation]["inter"][thread_count][0].append(comp_time)
            d[implementation]["inter"][thread_count][1].append(comm_time)  
    
    # Average all results for each tread count
    for implementation in d:
        for var in d[implementation]:
            for tc in d[implementation][var]:
                d[implementation][var][tc][0] = np.mean(d[implementation][var][tc][0]) # intra
                d[implementation][var][tc][1] = np.mean(d[implementation][var][tc][1]) # inter


    # Thread-counts should be equal between intra and inter node list for best comparison
    # This is limited by intra-node thread counts
    for implementation in d:
        intra_threads: List[int] = list(d[implementation]["intra"].keys())
        inter_threads: List[int] = list(d[implementation]["inter"].keys())
        if intra_threads != inter_threads:
            print(f"WARNING: Thread counts between inter and intra node lists not equal for {implementation}!", file=sys.stderr)

    # TODO: check if sorted on thread counts (should be by default, but still)
            
    if average:
        d = average_timings(d)

    # Set correct typing
    for implementation in d:
        for var in d[implementation]:
            comp_times: List[float] = []
            comm_times: List[float] = []
            threads: List[int] = []
            for tc in d[implementation][var]:
                threads.append(tc)
                comp_times.append(d[implementation][var][tc][0])
                comm_times.append(d[implementation][var][tc][1])
            d[implementation][var] = [comp_times, comm_times, threads]

    return d

def get_scaling_factors(data: Dict[str, Dict[str, Dict[int, List[float]]]]) -> Dict[str, Dict[str, Dict[int, List[float]]]]:
    """ Replaces the timing lists from data with scaling factor lists. 
        Scaling factor is determined by dividing comp/comm time at index i by their base value (index 0).
        The lowest thread count receives scaling factor 1.
        It is recommended to start with 2 threads.
        Otherwise, results will be distorted by (near-)zero communication time with 1 thread.""" 
    for implementation in data:
        for v in data[implementation]:
            # v = (intra || inter)
            comp_time: List[float] = data[implementation][v][0]
            comm_time: List[float] = data[implementation][v][1]
            threads: List[int] = data[implementation][v][2]
            if len(threads) != len(comp_time) or len(threads) != len(comm_time):
                print(f"ERROR: length mismatch between threads and comp/comm time lists!", file=sys.stderr)
                return data
            
            # Set base values for comp/comm time
            base_comp: float = comp_time[0]
            base_comm: float = comm_time[0]
            
            # Create scaling factor lists
            comp_factor: List[float] = [1]
            comm_factor: List[float] = [1]
            for i in range(len(threads)-1):
                index: int = i+1
                comp_factor.append(comp_time[index]/base_comp)
                comm_factor.append(comm_time[index]/base_comm)
            
            # Substitute timing lists with scaling factor lists
            data[implementation][v][0] = comp_factor
            data[implementation][v][1] = comm_factor
    
    # pprint(data)

    return data

def scaling_factor_plot(data: Dict[str, List[List[any]]], implementation: str, options: Dict[str, any]) -> None:
    """ Creates scaling factor plot comparing the scaling factors of communication and computation time. """
    datapath: str = options["datapath"]
    ex_num: int = options["ex-num"]
    cluster: str = options["compute-cluster"]
    date: str = options["file-date"]

    ### Intra node plot
    d = data["intra"]
    x: List[int] = d[2] # thread count list
    y_0: List[float] = d[0] # comp time scaling factors
    y_1: List[float] = d[1] # comm time scaling factors
    plt.plot(x, y_0, label="computation time")
    plt.plot(x, y_1, label="communication time")

    # plot parameters
    plt.xlabel("Thread count")
    plt.ylabel("Scaling factor")
    plt.legend(loc="upper right", ncol=1)
    # plt.title("Intra-node computation vs communication scaling")
    plt.savefig(f"{datapath}/ex{ex_num}_{cluster}_{implementation}_intra_{date}.png")
    plt.close()

    ### Inter node plot
    d = data["inter"]
    x: List[int] = d[2] # thread count list
    y_0: List[float] = d[0] # comp time scaling factors
    y_1: List[float] = d[1] # comm time scaling factors
    plt.plot(x, y_0, label="computation time")
    plt.plot(x, y_1, label="communication time")

    # plot parameters
    plt.xlabel("Thread count")
    plt.ylabel("Scaling factor")
    plt.legend(loc="upper right", ncol=1)
    # plt.title("Inter-node computation vs communication scaling")
    plt.savefig(f"{datapath}/ex{ex_num}_{cluster}_{implementation}_inter_{date}.png")
    plt.close()

    ### Intra vs Inter node plot
    d_0 = data["intra"]
    d_1 = data["inter"]
    x: List[int] = d_0[2] # thread count list
    y_0: List[float] = d_0[1] # comm time scaling factors
    y_1: List[float] = d_1[1] # comm time scaling factors
    plt.plot(x, y_0, label="intra-node comm. time")
    plt.plot(x, y_1, label="inter-node comm. time")

    # plot parameters
    plt.xlabel("Thread count")
    plt.ylabel("Scaling factor")
    plt.legend(loc="upper right", ncol=1)
    # plt.title("Intra-node vs inter-node communication scaling")
    plt.savefig(f"{datapath}/ex{ex_num}_{cluster}_{implementation}_comm_intra_v_inter_{date}.png")
    plt.close()

def create_plots(data: Dict[str, Dict[str, List[List[any]]]], options: Dict[str, any]) -> None:
    """ Creates scaling factor plot for each implementation in data. """
    for implementation in data:
        plot_data = data[implementation]
        scaling_factor_plot(plot_data, implementation, options)

def main():
    # Argument parsing
    parser: ArgumentParser = ArgumentParser()
    common.add_parameters(parser)
    parser.add_argument('--average', action='store_true')
    args = parser.parse_args()
    options: Dict[str, any] = dict(vars(args))

    prefix = f"EX{options['ex-num']}-{options['compute-cluster']}-TIMING-RESULTS-"
    file: str = options["datapath"] + "/" + prefix + options["file-date"] + options["file-extension"]
    if not Path(file).is_file():
        print(f"Error: {file} is not a file.", file=sys.stderr)
        return 1
    
    # Preprocessing
    df: DataFrame = pd.read_csv(file, delim_whitespace=True, names=common.timing_headers)
    average: bool = True if args.average else False
    data: Dict[any, any] = preprocess(df, average)
    data = get_scaling_factors(data)

    # Plot data
    create_plots(data, options)


if __name__ == "__main__":
    sys.exit(main())