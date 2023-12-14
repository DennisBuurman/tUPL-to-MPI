#!/usr/bin/env python3

# 
# Script to visualize results produced by process-results.py.
# Creates a graph for each experiment.
# 
# Author: Dennis Buurman, Leiden University

import sys
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import numpy as np

from typing import Dict, Tuple, List

filename = "EX{experiment_number}-{cluster}-RESULTS.txt"
path = "tupl-kmeans-39f1073/support/exp-{experiment_number}-{cluster}/"

# Column headers of results file entries 
headers = ["Implementation", "Input Size", "Clusters", "Dimension", "Nodes", "Tasks Per Node", "Worst Calc. Time", "Worst Exec. Time", "Iterations"]

# Implementation number as compared to Annes work
names = {
    "own": "Implementation 1",
    "own_inc": "Implementation 2",
    "own_loc": "Implementation 3",
    "own_inc_loc": "Implementation 4"
}

def process_data_exp1(df: DataFrame) -> Dict:
    d: Dict = {}
    df.reset_index()

    # Sort data on input size, implementation, worst calc. time
    for index, row in df.iterrows():
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
           implementations[i]: float = sum(times) / len(times)
    
    return d

def serialize_exp1(data: Dict) -> Tuple[List, Dict]:
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

def create_plot_exp1(sizes: List[int], times: Dict[str, float], name="") -> None:
    x = np.arange(len(sizes))
    width = 0.2 # bar width
    multiplier = 0

    fig, ax = plt.subplots(layout="constrained")

    for implementation, time in times.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, time, width, label=implementation)
        # ax.bar_label(rects, padding=3)
        multiplier += 1
    
    ax.set_ylabel("Y")
    ax.set_title("Input size variation")
    ax.set_xticks(x + width, sizes)
    ax.legend(loc="upper left", ncols=1)
    # ax.set_ylim(0, 3)

    plt.savefig(f"experiment_1_{name}.png")
    # plt.show()

# Varying input sizes
def exp1() -> None:
    # DAS-5 only
    cluster = "DAS5"
    file = path.format(experiment_number=1, cluster=cluster) + filename.format(experiment_number=1, cluster="DAS5")
    df = pd.read_csv(file, delim_whitespace=True, names=headers)
    data = process_data_exp1(df)
    sizes, times = serialize_exp1(data)
    create_plot_exp1(sizes, times, name=cluster)
    

    # DAS-6 only
    cluster = "DAS6"
    # TODO: repeat for DAS-6

    # Combined (original DAS-4?)
    # TODO: normalize DAS-5 and DAS-6 results (maybe add DAS-4 results)
    # TODO: compare normalized results in one plot

# Varying thread counts
def exp2() -> None:
    # TODO
    pass

# Varying dimensions
def exp3() -> None:
    pass

# Varying cluster sizes
def exp4() -> None:
    pass

def main():
    exp1()

    return 0

if __name__ == "__main__":
    sys.exit(main())