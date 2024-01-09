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
from datetime import datetime

from typing import Dict, Tuple, List

# Date string
date = datetime.today().strftime("%d-%m-%Y")

# Column headers of results file entries 
headers = ["Implementation", "Input Size", "Clusters", "Dimension", "Nodes", "Tasks Per Node", "Worst Calc. Time", "Worst Exec. Time", "Iterations"]

# Implementation number as compared to Annes work
names = {
    "own": "Implementation 1",
    "own_inc": "Implementation 2",
    "own_loc": "Implementation 3",
    "own_inc_loc": "Implementation 4"
}

def process(df: DataFrame) -> Dict:
    """ Processes data from a result file into a dictionary.
        Preps data for general graph. """
    d: Dict = {}
    df.reset_index()

    # Sort data on input size, implementation, worst calc. time
    for _, row in df.iterrows():
        implementation: str = row["Implementation"]
        input_size: str = row["Input Size"]
        time: float = row["Worst Calc. Time"]
        
        if input_size in d:
            time_per_implementation: Dict = d[input_size]
            if implementation in time_per_implementation:
                time_per_implementation[implementation].append(time)
            else:
                time_per_implementation[implementation] = [time]
        else:
            d[input_size]: Dict = {
                implementation: [
                    time
                ]
            }

    # Get average time per size per implementation
    for size in d:
        implementations: Dict = d[size]
        for i in implementations:
           times: List[float] = implementations[i]
           implementations[i]: float = np.mean(times)
    
    return d

def serialize(data: Dict) -> Tuple[List, Dict]:
    """ Serializes the data for matplotlib barplot input. """
    sizes: List[int] = []
    times: Dict = {}

    for size in data:
        sizes.append(size)
        for implementation in data[size]:
            t: float = data[size][implementation]
            n: str = names[implementation] # Anne's naming convention
            if n in times:
                times[n].append(t)
            else:
                times[n] = [t]

    return sizes, times

def create_plot(sizes: List[int], times: Dict[str, float], options) -> None:
    """ Creates the default plot for experiment 1: input size variation. """
    name = options["compute-cluster"] + "_" + date

    x = np.arange(len(sizes))
    width = 0.2 # bar width
    multiplier = 0

    fig, ax = plt.subplots(layout="constrained")

    for implementation, time in times.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, time, width, label=implementation)
        # ax.bar_label(rects, padding=3)
        multiplier += 1
    
    ax.set_ylabel("Calculation time (s)")
    ax.set_xlabel("Input size (2^x)")
    ax.set_title(f"Input size variation: {name}")
    ax.set_xticks(x + width, sizes)
    ax.legend(loc="upper left", ncols=1)
    # ax.set_ylim(0, 3)

    plt.savefig(f"exp_1_{name}.png")
    plt.close(fig)
    # plt.show()

def process_reverse(df: DataFrame) -> Dict:
    """ Processes data from a result file into a dictionary.
        Preps data for boxplot graphs. 
        Reverses implementation and input size compared to process().
        Removes the average time, as all results are needed for the boxplots. """
    d: Dict = {}
    df.reset_index()

    # Sort data on implementation, input size, worst calc. time
    for index, row in df.iterrows():
        implementation: str = row["Implementation"]
        input_size: str = row["Input Size"]
        time: float = row["Worst Calc. Time"]
        
        if implementation in d:
            time_per_size: Dict = d[implementation]
            if input_size in time_per_size:
                time_per_size[input_size].append(time)
            else:
                time_per_size[input_size] = [time]
        else:
            d[implementation]: Dict = {
                input_size: [
                    time
                ]
            }
    
    return d

def create_boxplots(data: Dict, options) -> None:
    """ Creates boxplots for each implementation to visualize performance variability. """
    cluster = options["compute-cluster"]

    for implementation in data:
        fig, ax = plt.subplots(layout="constrained")
        box_plot_data = []
        labels = []
        for size in data[implementation]:
            times = data[implementation][size]
            box_plot_data.append(times)
            labels.append(size)
        ax.set_title(f"{names[implementation]} {cluster} exp. 1 runtime distribution {size} ({date})")
        ax.set_ylabel("Calculation time (s)")
        ax.set_xlabel("Input size (2^x)")
        ax.boxplot(box_plot_data, patch_artist=True, labels=labels)
        plt.savefig(f"exp_1_bp_{implementation}_{cluster}_{date}")
        plt.close(fig)

def calc_ci(values, z=1.96) -> Tuple[float, float]:
    """ Calculates the 95% CI for given input. """
    mean = np.mean(values)
    stdev = np.std(values)
    ci = z * stdev / np.sqrt(len(values))
    return mean, ci

def ci_list(values) -> Tuple[List[float], List[float]]:
    """ Returns a list of 95% CI values. """
    mean_list = []
    ci_list = []
    for i in range(len(values)):
        mean, ci = calc_ci(values[0:i+1])
        mean_list.append(mean)
        ci_list.append(ci)
    return mean_list, ci_list

def create_confidence_interval(data: Dict, options) -> None:
    """ Creates a 95% CI plot over the runs of each implementation. """
    cluster = options["compute-cluster"]

    size: int = 28 # denotes size to create CI over
    for implementation in data:
        fig, ax = plt.subplots(layout="constrained")
        times: List[float] = data[implementation][size]
        x = np.arange(1, len(times)+1)
        mean_values, ci_values = ci_list(times)
        
        ax.set_title(f"{names[implementation]} {cluster} 95% CI on size {size} ({date})")
        ax.set_ylabel("Mean calculation time (s)")
        ax.set_xlabel("Iteration")
        ax.plot(x, mean_values)
        ax.fill_between(x, np.subtract(mean_values, ci_values), np.add(mean_values, ci_values), color='b', alpha=.1)

        plt.savefig(f"exp_1_ci_{implementation}_{cluster}_{size}_{date}")
        plt.close(fig)

def main():
    # Argument parsing
    parser = ArgumentParser()
    parser.add_argument("--compute-cluster", dest="compute-cluster", type=str, default="DAS5",
                            help="Compute cluster of the results")
    parser.add_argument("--file-preamble", dest="file-preamble", type=str, default="EX1-DAS5-RESULTS-",
                            help="Preamble of results file")
    parser.add_argument("--file-date", dest="file-date", type=str, default=date,
                            help="Date of the results file")
    parser.add_argument("--file-extension", dest="file-extension", type=str, default=".txt",
                            help="File extension of results file")
    parser.add_argument("--datadir", dest="datadir", type=str, default="",
                            help="Directory of the results file; include '/' at the end if non-empty")
    args = parser.parse_args()
    options = dict(vars(args))

    # Preprocessing
    file: str = options["datadir"] + options["file-preamble"] + options["file-date"] + options["file-extension"]
    if not Path(file).is_file():
        print(f"Error: {file} is not a file.", file=sys.stderr)
        return 1
    df: DataFrame = pd.read_csv(file, delim_whitespace=True, names=headers)
    data: Dict = process(df)
    sizes, times = serialize(data)

    # Main plot
    create_plot(sizes, times, options)
    # Box plots
    data = process_reverse(df)
    create_boxplots(data, options)
    # Confidence intervals
    create_confidence_interval(data, options)

    return 0

if __name__ == "__main__":
    sys.exit(main())