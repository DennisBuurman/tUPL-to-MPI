#!/usr/bin/env python3

# 
# Script to visualize results produced by process-results.py.
# Created to produce graphs for input size variation experiment;
# Also denoted as experiment 1.
# 
# Author: Dennis Buurman, Leiden University

import sys
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import numpy as np

from typing import Dict, Tuple, List

date = "20-12-2023" # DD-MM-YYYY
cluster = "DAS5" # DAS5 || DAS6
filename = f"EX2-{cluster}-RESULTS.txt"
path = f"{cluster}/EX2/{date}/"

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
        thread_count: int = row["Nodes"] * row["Tasks Per Node"]
        time: float = row["Worst Calc. Time"]
        
        if thread_count in d:
            time_per_implementation: Dict = d[thread_count]
            if implementation in time_per_implementation:
                time_per_implementation[implementation].append(time)
            else:
                time_per_implementation[implementation] = [time]
        else:
            d[thread_count]: Dict = {
                implementation: [
                    time
                ]
            }

    # Get average time per size per implementation
    for thread_count in d:
        implementations: Dict = d[thread_count]
        for i in implementations:
           times: List[float] = implementations[i]
           implementations[i]: float = np.mean(times)
    
    return d

def serialize(data: Dict) -> Tuple[List, Dict]:
    """ Serializes the data for matplotlib barplot input. """
    thread_count_list: List[int] = []
    times: Dict = {}

    for thread_count in data:
        if thread_count > 500:
            continue
        thread_count_list.append(thread_count)
        for implementation in data[thread_count]:
            t: float = data[thread_count][implementation]
            n: str = names[implementation] # Anne's naming convention
            if n in times:
                times[n].append(t)
            else:
                times[n] = [t]

    return thread_count_list, times

def create_plot(thread_count: List[int], times: Dict[str, float], name="") -> None:
    """ Creates the default plot for experiment 2: thread count variation. """
    x = np.arange(len(thread_count))

    fig, ax = plt.subplots(layout="constrained")
    markers = "x+*2"
    count = 0

    for implementation, time in times.items():
        plt.plot(thread_count, time, marker=markers[count], label=implementation)
        count += 1
    
    # plt.yscale("log")

    ax.set_ylabel("Calculation time (s)")
    ax.set_xlabel("Thread count")
    ax.set_title(f"Thread count variation: {name}")
    ax.set_xticks(thread_count)
    ax.legend(loc="upper right", ncols=1)
    # ax.set_ylim(0, 3)

    plt.savefig(path+f"exp_2_{name}.png")
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
        thread_count: int = row["Nodes"] * row["Tasks Per Node"]
        time: float = row["Worst Calc. Time"]
        
        if implementation in d:
            time_per_size: Dict = d[implementation]
            if thread_count in time_per_size:
                time_per_size[thread_count].append(time)
            else:
                time_per_size[thread_count] = [time]
        else:
            d[implementation]: Dict = {
                thread_count: [
                    time
                ]
            }
    
    return d

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

def create_confidence_interval(data: Dict) -> None:
    """ Creates a 95% CI plot over the runs of each implementation. """
    thread_count: int = 320 # denote thread count to create CI over
    for implementation in data:
        fig, ax = plt.subplots(layout="constrained")
        times: List[float] = data[implementation][thread_count]
        x = np.arange(1, len(times)+1)
        mean_values, ci_values = ci_list(times)
        
        ax.set_title(f"{names[implementation]} {cluster} {thread_count} threads 95% CI ({date})")
        ax.set_ylabel("Mean calculation time (s)")
        ax.set_xlabel("Iteration")
        ax.plot(x, mean_values)
        ax.fill_between(x, np.subtract(mean_values, ci_values), np.add(mean_values, ci_values), color='b', alpha=.1)

        plt.savefig(path+f"exp_2_ci_{implementation}_{cluster}_tc{thread_count} ({date})")
        plt.close(fig)

def main():
    file: str = path + filename.format(experiment_number=1, cluster=cluster)
    df: DataFrame = pd.read_csv(file, delim_whitespace=True, names=headers)
    data: Dict = process(df)
    thread_count, times = serialize(data)
    create_plot(thread_count, times, name=f"{cluster} ({date})")

    data = process_reverse(df)
    create_confidence_interval(data)

if __name__ == "__main__":
    sys.exit(main())