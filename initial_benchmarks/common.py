#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
import pandas as pd
from pandas import DataFrame
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

from typing import Dict, Tuple, List

# Column headers of results file entries 
headers: List[str] = [
    "Implementation", 
    "Input Size", 
    "Clusters", 
    "Dimension", 
    "Nodes", 
    "Tasks Per Node", 
    "Worst Calc. Time", 
    "Worst Exec. Time", 
    "Iterations"
]

# Implementation number as compared to Annes work
names: Dict[str, str] = {
    "own": "Implementation 1",
    "own_inc": "Implementation 2",
    "own_loc": "Implementation 3",
    "own_inc_loc": "Implementation 4",
    "own_m": "Implementation 5",
    "own_m_loc": "Implementation 6",
    "own_im": "Implementation 7",
    "own_im_loc": "Implementation 8",
    "own_values_only": "Implementation 9",
    "own_values_only_loc": "Implementation 10"
}

def add_parameters(parser: ArgumentParser) -> None:
    parser.add_argument("--datapath", dest="datapath", type=str, default="results",
                            help="Location of the results file to process")
    parser.add_argument("--compute-cluster", dest="compute-cluster", type=str, default="DAS5",
                            help="Compute cluster of the results")
    parser.add_argument("--file-date", dest="file-date", type=str, default=datetime.today().strftime("%d-%m-%Y"),
                            help="Date of the results file")
    parser.add_argument("--file-extension", dest="file-extension", type=str, default=".txt",
                            help="File extension of results file")
    parser.add_argument("--ex-num", dest="ex-num", type=int, default=0,
                            help="experiment number used for file naming in experiment scripts")

def process(df: DataFrame, variable: str) -> Dict[str, Dict[str, List[float]]]:
    """ Processes data from a result file into a dictionary.
        Preps data for general graph. 
        variable should be in the headers global list."""
    d: Dict[str, Dict[str, List[float]]] = {}
    df.reset_index()

    if variable not in list(df.columns):
        print(f"ERROR: variable '{variable}' not in dataframe!", file=sys.stderr)
        return {}
    
    # Sort data on variable, implementation, worst calc. time
    for _, row in df.iterrows():
        implementation: str = row["Implementation"]
        var: str = row[variable]
        time: float = row["Worst Calc. Time"]
        if var in d:
            time_per_implementation: Dict = d[var]
            if implementation in time_per_implementation:
                time_per_implementation[implementation].append(time)
            else:
                time_per_implementation[implementation] = [time]
        else:
            d[var] = {
                implementation: [
                    time
                ]
            }
    return d

def process_reverse(df: DataFrame, variable: str) -> Dict:
    """ Processes data from a result file into a dictionary.
        Preps data for boxplot and CI graphs. 
        Reverses implementation and input size compared to process(). """
    d: Dict[str, Dict[str, List[float]]] = {}
    df.reset_index()

    if variable not in list(df.columns):
        print(f"ERROR: variable '{variable}' not in dataframe!", file=sys.stderr)
        return {}

    # Sort data on implementation, variable, worst calc. time
    for _, row in df.iterrows():
        implementation: str = row["Implementation"]
        var: str = row[variable]
        time: float = row["Worst Calc. Time"]
        
        if implementation in d:
            time_per_var: Dict[str, Dict[str, List[float]]] = d[implementation]
            if var in time_per_var:
                time_per_var[var].append(time)
            else:
                time_per_var[var] = [time]
        else:
            d[implementation] = {
                var: [
                    time
                ]
            }
    
    return d

def average_calc_times(d: Dict[str, Dict[str, List[float]]]) -> Dict[str, Dict[str, float]]:
    """ Calculates the mean of the calculation times present in a nested dictionary created by process(). 
        List of float values is replaced by its mean. """
    for x in d:
        implementations: Dict[str, List[float]] = d[x]
        for i in implementations:
           times: List[float] = implementations[i]
           implementations[i] = np.mean(times)
    
    return d

def serialize(data: Dict[str, Dict[str, float]]) -> Tuple[List[int], Dict]:
    """ Serializes the data for matplotlib barplot input. """
    var_list: List[int] = [] # list of variable values
    times: Dict[str, List[float]] = {}

    for var in data:
        var_list.append(var)
        for implementation in data[var]:
            t: float = data[var][implementation]
            n: str = names[implementation] # Anne's naming convention
            if n in times:
                times[n].append(t)
            else:
                times[n] = [t]

    return var_list, times

def calc_ci(values: List[float], z: float = 1.96) -> Tuple[float, float]:
    """ Calculates the 95% CI for given input. """
    mean: float = np.mean(values)
    stdev: float = np.std(values)
    ci = z * stdev / np.sqrt(len(values))
    return mean, ci

def ci_list(values) -> Tuple[List[float], List[float]]:
    """ Returns a list of 95% CI values. """
    mean_list: List[float] = []
    ci_list: List[float] = []
    for i in range(len(values)):
        mean, ci = calc_ci(values[0:i+1])
        mean_list.append(mean)
        ci_list.append(ci)
    return mean_list, ci_list

def create_confidence_interval(data: Dict[str, Dict[str, List[float]]], options: Dict[str, any], variable: int, ex_num: int) -> None:
    """ Creates a 95% CI plot over the runs of each implementation. 
        'variable' denotes the variable value, of the var used in process functions, to create CI over. """
    cluster = options["compute-cluster"]
    datapath = options["datapath"]
    date = options["file-date"]

    for implementation in data:
        fig, ax = plt.subplots()
        fig.tight_layout()
        times: List[float] = data[implementation][variable]
        x = np.arange(1, len(times)+1)
        mean_values, ci_values = ci_list(times)
        
        ax.set_title(f"{names[implementation]} {cluster} 95% CI on size 2^{variable}")
        ax.set_ylabel("Mean calculation time (s)")
        ax.set_xlabel("Iteration")
        ax.plot(x, mean_values)
        ax.fill_between(x, np.subtract(mean_values, ci_values), np.add(mean_values, ci_values), color='b', alpha=.1)

        plt.savefig(f"{datapath}/ex{ex_num}_{cluster}_ci_{implementation}_size{variable}_{date}")
        plt.close(fig)
