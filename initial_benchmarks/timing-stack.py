#!/usr/bin/env python3

# 
# Script to visualize results produced by process-timing.py.
# Created to produce stacked bar plot graphs for time spent in components;
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

# Variables used for stacked bar plot
variables = [
    "Worst Init. Time", 
    "Worst Reassign Time",
    "Worst Comm. Time", 
    "Worst Recalc. Time",
]

def average_data(data: Dict[str, Dict[str, Dict[str, List[float]]]]) -> Dict[str, Dict[str, Dict[str, float]]]:
    """ Averages the timing variables present in the list of floats from data. """
    for size in data:
        for implementation in data[size]:
            for timing_variable in data[size][implementation]:
                l: List[float] = data[size][implementation][timing_variable]
                data[size][implementation][timing_variable] = np.mean(l)
    
    return data

def preprocess(df: DataFrame) -> Dict[str, Dict[str, Dict[str, List[float]]]]:
    """ Preprocessing for stacked bar plots. 
        Creates """
    d: Dict[str, Dict[str, Dict[str, List[float]]]] = dict()
    df.reset_index()

    for v in variables:
        if v not in list(df.columns):
            print(f"ERROR: variable '{v}' not in dataframe!", file=sys.stderr)
            return d
    
    # sort on size, implementation, timing variables
    # {size: {implementation: {time: []}}}
    for _, row in df.iterrows():
        size: int = row["Input Size"]
        implmentation: str = row["Implementation"]
        if size not in d:
            d[size] = {}
        if implmentation not in d[size]:
            d[size][implmentation] = {}
        for v in variables:
            v_time: float = row[v]
            if v not in d[size][implmentation]:
                d[size][implmentation][v] = list()
            d[size][implmentation][v].append(v_time)

    return d

def stacked_plot(data: Dict[str, Dict[str, float]], size: int):
    x: List[str] = list() # [implementations]
    Y: List[List[float]] = list() # list of times per implementation
    legend: List[str] = list() # list of names per timing variable
    
    # NOTE: Y needs to be rotated and reversed before a stacked plot can be created!
    for implementation in data:
        x.append(implementation)
        y: List[float] = list()
        for v in data[implementation]:
            y.append(data[implementation][v])
            if v not in legend:
                legend.append(v)
        Y.append(y)

    # Rotate and inverse Y to get 'stackable' lists
    Y = list(zip(*Y))[::-1][::-1]

    # Create plot
    for index in range(len(Y)):
        bottom: np.array = np.array([0.0] * len(Y[index]))
        if (index > 0):
            for i in range(index, 0, -1):
                bottom += np.array(Y[i-1])
        plt.bar(x, Y[index], bottom=bottom)

    # Add plot parameters
    plt.xlabel("Implementation")
    plt.ylabel("Time (s)")
    plt.legend(legend)
    plt.title("Average Iteration Time Distribution")
    # plt.savefig(f"{datapath}/ex{ex_num}_{cluster}_ci_{implementation}_size{variable}_{date}")
    plt.savefig(f"testing_s{size}")
    plt.close()
            

def create_plots(data: Dict[str, Dict[str, Dict[str, float]]]):
    """ Create plot for each size in data argument. """
    for size in data:
        plot_data = data[size]
        stacked_plot(plot_data, size)

def main():
    # Argument parsing
    parser: ArgumentParser = ArgumentParser()
    common.add_parameters(parser)
    args = parser.parse_args()
    options: Dict[str, any] = dict(vars(args))

    prefix = f"EX{options['ex-num']}-{options['compute-cluster']}-TIMING-RESULTS-"
    file: str = options["datapath"] + "/" + prefix + options["file-date"] + options["file-extension"]
    if not Path(file).is_file():
        print(f"Error: {file} is not a file.", file=sys.stderr)
        return 1
    
    # Preprocessing
    df: DataFrame = pd.read_csv(file, delim_whitespace=True, names=common.timing_headers)
    data: Dict[str, Dict[str, Dict[str, List[float]]]] = preprocess(df)
    data = average_data(data)

    # Plot data
    create_plots(data)


if __name__ == "__main__":
    sys.exit(main())