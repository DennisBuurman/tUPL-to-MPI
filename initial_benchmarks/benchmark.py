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

# Date string
date = datetime.today().strftime("%d-%m-%Y")

# Experiment 1 parameters
ex1_config = {
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
ex2_config = {
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

execs = ["own", "own_loc", "own_inc", "own_inc_loc"]

def run_experiment(config, options, ex_num):
    # Configuration variables
    variant = config['variant']
    size = config['size']
    clusters = " ".join([str(x) for x in config['clusters']])
    dimension = " ".join([str(x) for x in config['dimension']])
    seed = config['seed']
    nodes = config['nodes']
    tasks = config['ntasks-per-node']
    repeat = config['repeat']

    # Argument variables
    output_dir = options["outputdir"]
    compute_cluster = options["compute-cluster"]
    datapath = options["datapath"]
    scriptpath = options["scriptpath"]

    print("> Running experiment 1 ...")

    # Check if any .out files still in directory
    old_results = len(glob.glob1(".","*.out"))
    if old_results > 0:
        print(f"WARNING: {old_results} '*.out' files already present in directory!")

    # Prepare commands and execute script
    file = scriptpath + "/submit-exp.py"
    f = Path(file)
    if not f.is_file():
        print(f"ERROR: {file} is not a file.", file=sys.stderr)
        return 1

    if len(nodes) != len(tasks) or len(nodes) < 1:
        print(f"ERROR: nodes and tasks should be a list of list of int containing nodes * tasks configurations constructed using the index of the outer list. \
              Example: [[1], [8]] [[2, 4], [3]] refers to 1*2, 1*4, and 8*3")
        return 1

    # Save amount of node*task configs for progress checking
    config_counter = 0

    # Create commands for each nodes * tasks config
    for i in range(len(nodes)):
        n_list = nodes[i]
        t_list = tasks[i]
        node_str = ""
        for n in n_list:
            node_str += str(n) + " "
        tasks_str = ""
        for t in t_list:
            tasks_str += str(t) + " "
        command = f"./{file} --datapath {datapath} --variant {variant} --size {size} --clusters {clusters} --dimension {dimension} --seed {seed} --nodes {node_str} --ntasks-per-node {tasks_str} --repeat {repeat}"
        print(f"> Running commmand: {command}")
        subprocess.run(command.split())
        config_counter += len(n_list) * len(t_list)

    # Check progress of experiment
    filecount = config_counter * len(variant.split()) * len(size.split()) * len(clusters.split()) * len(dimension.split())
    current = 0
    new = len(glob.glob1(".","*.out"))
    pbar = tqdm(total=filecount)
    while (new - old_results < filecount):
        current = new
        time.sleep(1)
        new = len(glob.glob1(".","*.out"))
        pbar.update(new - current)
    pbar.close()
        
    # s_sleep = 300 # ~180s read time and 10*6s calc. time for size 28
    for i in tqdm(range(0, s_sleep), total = s_sleep, desc ="Waiting for jobs to finish"):
        time.sleep(1)

    # Process results
    print("> Processing available results ...")
    if not Path(output_dir).is_dir():
        Path(output_dir).mkdir(parents=True)
    filename = f"{output_dir}/EX1-{compute_cluster}-RESULTS-{date}.txt"
    result_file = open(filename, "w")
    file = scriptpath + "/process-results.py"
    f = Path(file)
    if not f.is_file():
        print(f"Error: {file} is not a file.", file=sys.stderr)
        return 1
    subprocess.run([f"./{file}", "."], stdout=result_file)
    print(f"> Running command: ./{file} .")
    print(f"> Results saved in {filename}")
    result_file.close()

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
    args = parser.parse_args()
    options = dict(vars(args))

    # Run selected experiment
    experiment = options["experiment"]
    del options["experiment"]
    if experiment == 1:
        return(run_experiment(ex1_config, options, 1))
    elif experiment == 2:
        return(run_experiment(ex2_config[options["compute-cluster"]], options, 2))

if __name__ == "__main__":
    sys.exit(main())