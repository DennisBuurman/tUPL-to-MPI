#!/usr/bin/env python3

# 
# Script to run the benchmark of experiment 1 and 2 of Anne Hommelbergs work.
# Results will be processed into .txt files, which can be visualized using
# the corresponding visualize.py file.
# 
# Author: Dennis Buurman, Leiden University

from pathlib import Path
import os
import sys
from argparse import ArgumentParser
import subprocess
from datetime import datetime

from typing import Dict

# Date string
date = datetime.today().strftime("%d-%m-%Y")

# Path to submit-exp.py script
supportdir = "../tupl-kmeans-39f1073/support/"

# Default settings
defaults: Dict = {"outputdir": "results",
                  "datapath": "/var/scratch/dbuurman/kmeans",
                  "seed": 971,
                  "compute-cluster": "DAS5",
                  "experiment": 1,
                  "variants": ["own", "own_inc", "own_loc", "own_inc_loc"],
                  "sizes": [20, 21, 22, 23, 24, 25, 26, 27, 28],
                  "clusters": [4],
                  "dimensions": [4],
                  "nodes": [8],
                  "ntasks-per-node": [8],
                  "repeats": 10}

# Implementation variants
execs = [
    "own",
    "own_loc",
    "own_inc",
    "own_inc_loc"
]

# Experiment variants
experiments = [1, 2]

def experiment_1(options) -> int:
    # Set options
    output_dir = options["outputdir"]
    compute_cluster = options["compute-cluster"]
    datapath = options["datapath"]
    variant = options["variant"] if type(options["variant"]) == str else " ".join(options["variant"])
    size = options["size"] if type(options["size"]) == str else " ".join([str(x) for x in options["size"]])
    clusters = options["clusters"] if type(options["clusters"]) == str else " ".join([str(x) for x in options["clusters"]])
    dimension = options["dimension"] if type(options["dimension"]) == str else " ".join([str(x) for x in options["dimension"]])
    seed = options["seed"]
    nodes = options["nodes"]
    tasks = options["ntasks-per-node"]
    repeats = options["repeats"]

    print("> Running experiment 1 ...")

    # Prepare command and execute script
    scriptpath = supportdir + "submit-exp.py"
    file = Path(scriptpath)
    if not file.is_file():
        print(f"Error: {scriptpath} is not a file.", file=sys.stderr)
        return 1
    command = f"./{scriptpath} --datapath {datapath} --variant {variant} --size {size} \
        --cluster {clusters} --dimension {dimension} --seed {seed} --nodes {nodes} --ntasks-per-node {tasks} --repeats {repeats}"
    print(f"--> Executing: {command}")
    subprocess.run(command.split())

    print("> Processing experiment 1 results ...")

    # Process results
    result_file = open(f"{output_dir}/EX1-{compute_cluster}-RESULTS-{date}.txt", "w")
    scriptpath = supportdir + "process-results.py"
    file = Path(scriptpath)
    if not file.is_file():
        print(f"Error: {scriptpath} is not a file.", file=sys.stderr)
        return 1
    subprocess.run([f"./{scriptpath}", "."], stdout=result_file)
    print(f"--> Executing: ./{scriptpath} .")
    result_file.close()

    print("> Visualizing experiment 1 results ...")

    # Visualize experiment results
    command = f"./ex1.py --compute-cluster {compute_cluster} --file-date {date} --datapath {output_dir}"
    print(f"--> Executing: {command}")
    subprocess.run(command.split())

    print("> Cleaning output files ...")

    # TODO

    print("> Done!")

    return 0

def experiment_2(options) -> int:
    print(options)

    # Set options
    output_dir = options["outputdir"]
    compute_cluster = options["compute-cluster"]
    datapath = options["datapath"]
    variant = options["variant"] if type(options["variant"]) == str else " ".join(options["variant"])
    size = options["size"] if type(options["size"]) == str else " ".join([str(x) for x in options["size"]])
    clusters = options["clusters"] if type(options["clusters"]) == str else " ".join([str(x) for x in options["clusters"]])
    dimension = options["dimension"] if type(options["dimension"]) == str else " ".join([str(x) for x in options["dimension"]])
    seed = options["seed"]
    nodes = options["nodes"]
    tasks = options["ntasks-per-node"]
    repeats = options["repeats"]

    print("> Running experiment 2 ...")

    # Prepare command and execute script
    scriptpath = supportdir + "submit-exp.py"
    file = Path(scriptpath)
    if not file.is_file():
        print(f"Error: {scriptpath} is not a file.", file=sys.stderr)
        return 1

    # List of tuples containing node and thread counts
    if compute_cluster == "DAS5":
        config = [([1], [2, 4, 8, 12, 16]),
                  ([1, 2], [12]),
                  ([1, 2, 3, 4, 5, 6, 7, 8, 12, 16, 20], [16])]
    elif compute_cluster == "DAS6":
        config = [([1], [2, 4, 8, 12, 16, 24]),
                  ([1, 2], [32]),
                  ([1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16], [48])]
    
    # Execute command for each node, thread tuple
    for tup in config:
        node_count = ""
        thread_count = ""
        for n in tup[0]:
            node_count += str(n) + " "
        for t in tup[1]:
            thread_count += str(t) + " "
        command = f"./{scriptpath} --datapath {datapath} --variant {variant} --size {size} \
                    --cluster {clusters} --dimension {dimension} --seed {seed} --nodes {node_count} --ntasks-per-node {thread_count} --repeats {repeats}"
        print(f"--> Executing: {command}")
        subprocess.run(command.split())


    print("> Processing experiment 2 results ...")

    # Process results
    result_file = open(f"{output_dir}/EX2-{compute_cluster}-RESULTS-{date}.txt", "w")
    scriptpath = supportdir + "process-results.py"
    file = Path(scriptpath)
    if not file.is_file():
        print(f"Error: {scriptpath} is not a file.", file=sys.stderr)
        return 1
    print(f"--> Executing: ./{scriptpath} .")
    subprocess.run([f"./{scriptpath}", "."], stdout=result_file)
    result_file.close()
    print(f"--> Result saved in 'EX2-{compute_cluster}-RESULTS-{date}.txt'")

    print("> Visualizing experiment 2 results ...")

    # Visualize experiment results
    command = f"./ex2.py --compute-cluster {compute_cluster} --file-date {date} --datapath {output_dir}"
    print(f"--> Executing: {command}")
    subprocess.run(command.split())

    print("> Cleaning output files ...")

    # TODO

    print("> Done!")

    return 0

def main() -> int:
    parser = ArgumentParser()

    # Parser arguments
    parser.add_argument("--outputdir", dest="outputdir", type=str, default=defaults["outputdir"],
                            help="Location of the resulting output")
    parser.add_argument("--datapath", dest="datapath", type=str, default=defaults["datapath"],
                            help="Location of the dataset")
    parser.add_argument("--seed", dest="seed", nargs="*", type=int, default=defaults["seed"],
                            help="Random seed the database was generated with")
    parser.add_argument("--compute-cluster", dest="compute-cluster", type=str, default=defaults["compute-cluster"],
                            help="Compute cluster in use")
    parser.add_argument("--experiment", dest="experiment", type=int, default=defaults["experiment"],
                            help="Experiment number to execute")
    parser.add_argument("--variant", dest="variant", nargs="*", default=" ".join(defaults["variants"]), choices=execs,
                        help="Variant of the algorithm to run")
    parser.add_argument("--size", dest="size", nargs="*", default=" ".join([str(x) for x in defaults["sizes"]]),
                            help="Size of the dataset")
    parser.add_argument("--clusters", dest="clusters", nargs="*", type=str, default=" ".join([str(x) for x in defaults["clusters"]]),
                            help="The number of clusters (k) to find")
    parser.add_argument("--dimension", dest="dimension", nargs="*", default=" ".join([str(x) for x in defaults["dimensions"]]),
                            help="The dimension of the data points")
    parser.add_argument("--nodes", dest="nodes", nargs="*", default=" ".join([str(x) for x in defaults["nodes"]]),
                        help="Number of nodes")
    parser.add_argument("--ntasks-per-node", dest="ntasks-per-node", nargs="*", default=" ".join([str(x) for x in defaults["ntasks-per-node"]]),
                        help="Number tasks to start on each node")
    parser.add_argument("--repeat", dest="repeat", type=int,
                        nargs=1, default=defaults["repeats"],
                        help="Number of times to repeat the experiment")

    # Prepare arguments
    args = parser.parse_args()
    options = dict(vars(args))
    experiment = options["experiment"]
    del options["experiment"]

    # Run experiment
    if experiment == 1:
        return(experiment_1(options))
    elif experiment == 2:
        return(experiment_2(options))

if __name__ == "__main__":
    sys.exit(main())