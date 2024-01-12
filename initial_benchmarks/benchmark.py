#!/usr/bin/env python3

# 
# Script to run the benchmark of experiment 1 and 2 of Anne Hommelbergs work.
# Results will be processed into .txt files, which can be visualized using
# the corresponding visualize.py file.
# 
# Author: Dennis Buurman, Leiden University
import sys
from argparse import ArgumentParser
import subprocess
from datetime import datetime
from pathlib import Path
import glob
from tqdm.auto import tqdm
import time

from typing import Dict, List

# Date string
date = datetime.today().strftime("%d-%m-%Y")

# Experiment 1 parameters
ex1_config: Dict[str] = {
    "seed": "971",
    "variant": "own_inc_loc",
    "size": "28",
    "clusters": [4],
    "dimension": [4],
    "nodes": [[8]],
    "ntasks-per-node": [[8]],
    "repeat": "10"
}

# Experiment 2 parameters
ex2_config: Dict[str] = {
    "DAS5": {
        "seed": "971",
        "variant": "own own_inc own_loc own_inc_loc",
        "size": "26",
        "clusters": [4],
        "dimension": [4],
        "nodes": [[1], [1, 2], [1, 2, 3, 4, 5, 6, 7, 8, 12, 16, 20]],
        "ntasks-per-node": [[2, 4, 8, 12, 16], [12], [16]],
        "repeat": "10"
    },
    "DAS6": {
        "seed": "971",
        "variant": "own own_inc own_loc own_inc_loc",
        "size": "26",
        "clusters": [4],
        "dimension": [4],
        "nodes": [[1], [1, 2], [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16]],
        "ntasks-per-node": [[2, 4, 8, 12, 16, 24], [32], [48]],
        "repeat": "10"
    }
}

execs: List[str] = ["own", "own_loc", "own_inc", "own_inc_loc"]

# Function to verify files
def exists(file: str) -> bool:
    f = Path(file)
    if not f.is_file():
        print(f"ERROR: {file} is not a file.", file=sys.stderr)
        return False
    return True

def submit_jobs(file: str, datapath: str, variant: str, size: str, clusters: str, dimension: str, seed: str, nodes: List[List[int]], tasks: List[List[int]], repeat: str) -> int:
    config_counter: int = 0 # count node*task configs
    for i in range(len(nodes)):
        n_list: List[int] = nodes[i]
        t_list: List[int] = tasks[i]
        node_str: str = ""
        for n in n_list:
            node_str += str(n) + " "
        tasks_str: str = ""
        for t in t_list:
            tasks_str += str(t) + " "
        command: str = f"./{file} --datapath {datapath} --variant {variant} --size {size} --clusters {clusters} --dimension {dimension} --seed {seed} --nodes {node_str} --ntasks-per-node {tasks_str} --repeat {repeat}"
        print(f"> Running commmand: {command}")
        subprocess.run(command.split())
        config_counter += len(n_list) * len(t_list)    
    return config_counter

def progress(filecount: int, old_results: int) -> None:
    current: int = 0
    new: int = len(glob.glob1(".","*.out"))
    pbar = tqdm(total=filecount)
    while (new - old_results < filecount):
        current = new
        time.sleep(1)
        new = len(glob.glob1(".","*.out"))
        pbar.update(new - current)
    pbar.close()

def sleep(seconds: int) -> None:
    for i in tqdm(range(0, seconds), total = seconds, desc ="Waiting for jobs to finish"):
        time.sleep(1)

def process_results(output_dir: str, compute_cluster: str, scriptpath: str) -> int:
    if not Path(output_dir).is_dir():
        Path(output_dir).mkdir(parents=True)
    filename = f"{output_dir}/EX1-{compute_cluster}-RESULTS-{date}.txt"
    result_file = open(filename, "w")
    file = scriptpath + "/process-results.py"
    if not exists(file):
        return 1
    subprocess.run([f"./{file}", "."], stdout=result_file)
    print(f"> Running command: ./{file} .")
    print(f"> Results saved in {filename}")
    result_file.close()
    return 0

def validate_results() -> None:
    # TODO
    pass

def report_results() -> None:
    # TODO
    pass

def run_experiment(config: Dict[str], options: Dict[str], ex_num: int) -> int:
    # Configuration variables
    variant: str = config['variant']
    size: str = config['size']
    clusters: str = " ".join([str(x) for x in config['clusters']])
    dimension: str = " ".join([str(x) for x in config['dimension']])
    seed: str = config['seed']
    nodes: List[List[int]] = config['nodes']
    tasks: List[List[int]] = config['ntasks-per-node']
    repeat: str = config['repeat']
    # Argument variables
    output_dir: str = options["outputdir"]
    compute_cluster: str = options["compute-cluster"]
    datapath: str = options["datapath"]
    scriptpath: str = options["scriptpath"]

    print(f"> Running experiment {ex_num} ...")

    # Check if any .out files still in directory
    old_results: int = len(glob.glob1(".","*.out"))
    if old_results > 0:
        print(f"WARNING: {old_results} '*.out' files already present in directory!")
    # Prepare commands and execute script
    file = scriptpath + "/submit-exp.py"
    if not exists(file):
        return 1
    # Check nodes and tasks lists
    if len(nodes) != len(tasks) or len(nodes) < 1:
        print(f"ERROR: nodes and tasks should be a list of list of int containing nodes * tasks configurations constructed using the index of the outer list. \
              Example: [[1], [8]] [[2, 4], [3]] refers to 1*2, 1*4, and 8*3")
        return 1
    # Create commands for each nodes * tasks config
    config_counter = submit_jobs(file, datapath, variant, size, clusters, dimension, seed, nodes, tasks, repeat)
    # Check progress of experiment
    filecount = config_counter * len(variant.split()) * len(size.split()) * len(clusters.split()) * len(dimension.split())
    progress(filecount, old_results)
    # Wait for results to finish
    sleep(300) # ~180s read time and 10*6s calc. time for size 28
    # Process results
    print("> Processing available results ...")
    res: int = process_results(output_dir, compute_cluster, scriptpath)
    if res != 0:
        return res
    # Create results report
    report_results()
    # Finish
    print("Done!")
    print(f"Visualize results by running:\n./ex{ex_num}.py --compute-cluster {compute_cluster} --file-date {date} --datapath {output_dir}")

    return 0

def main():
    parser = ArgumentParser()
    parser.add_argument("--compute-cluster", "-c", dest="compute-cluster", type=str, required=True, choices=["DAS5", "DAS6"],
                            help="Compute cluster in use")
    parser.add_argument("--experiment", "-e", dest="experiment", type=int, default=1, choices=[1, 2],
                            help="Experiment number to execute")
    parser.add_argument("--outputdir", "-o", dest="outputdir", type=str, default="results",
                            help="Location of the resulting output")
    parser.add_argument("--datapath", "-d", dest="datapath", type=str, default="/var/scratch/dbuurman/kmeans",
                            help="Location of the dataset")
    parser.add_argument("--scriptpath", "-s", dest="scriptpath", type=str, default="../tupl-kmeans-39f1073/support",
                            help="Location to directory of submit-exp.py script")
    # TODO: only validate and report existing results
    # TODO: clean *.out files argument
    args = parser.parse_args()
    options: Dict[str] = dict(vars(args))

    # Run selected experiment
    experiment: int = options["experiment"]
    del options["experiment"]
    if experiment == 1:
        return(run_experiment(ex1_config, options, 1))
    elif experiment == 2:
        return(run_experiment(ex2_config[options["compute-cluster"]], options, 2))

if __name__ == "__main__":
    sys.exit(main())