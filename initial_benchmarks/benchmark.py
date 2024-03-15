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
import os
import random
import operator

from typing import Dict, List, Tuple

# DAS account name
account_name: str = "dbuurman"

# Date string
date: str = datetime.today().strftime("%d-%m-%Y")

# Experiment 1 default parameters
ex1_config: Dict[str, any] = {
    "seed": "971",
    "variant": "own own_inc own_loc own_inc_loc",
    "size": "24 25 26 27 28",
    # "size": "24 25 26 27 28",
    # "size": "24",
    "clusters": [4],
    "dimension": [4],
    "nodes": [[8]],
    # "nodes": [[1]],
    "ntasks-per-node": [[8]],
    "repeat": "10"
}

# Experiment 2 default parameters
ex2_config: Dict[str, any] = {
    "DAS5": {
        "seed": "971",
        "variant": "own own_inc own_loc own_inc_loc",
        "size": "28",
        "clusters": [4],
        "dimension": [4],
        "nodes": [[1], [1, 2, 3, 4, 5, 6, 7, 8, 12, 16]],
        "ntasks-per-node": [[2, 4, 8, 12, 16, 24], [32]],
        "repeat": "10"
    },
    "DAS6": {
        "seed": "971",
        "variant": "own own_inc own_loc own_inc_loc",
        "size": "28",
        "clusters": [4],
        "dimension": [4],
        "nodes": [[1], [1, 2], [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16]],
        "ntasks-per-node": [[2, 4, 8, 12, 16, 24], [32], [48]],
        "repeat": "10"
    }
}

execs: List[str] = ["own", "own_inc", "own_loc", "own_inc_loc"]
debug = False

def exists(file: str) -> bool:
    """ Checks if given filename is a file in the current directory """
    f = Path(file)
    if not f.is_file():
        print(f"ERROR: {file} is not a file.", file=sys.stderr)
        return False
    return True

def split_arguments(arguments: List[str]) -> Dict[str, List[str]]:
    """ Splits a list of arguments into a dictionary with argument names as keys and argument values as values. """
    names: Dict[str, List[str]] = {}
    # Get arguments and values in dict
    for a in arguments:
        options = a.strip("\n").split()
        n = options[0] # argument name
        names[n] = options[1:] # argument options / values
    return names

def submit_jobs(file: str, datapath: str, variant: str, size: str, clusters: str, dimension: str, seed: str, nodes: List[List[int]], tasks: List[List[int]], repeat: str, compute_cluster: str, mean_set: int, init_seed: int) -> Tuple[int, List[str]]:
    """ Submits jobs to DAS5 or DAS6 cluster depending on configuration. Configurations must adhere ex1_config notation. """
    config_counter: int = 0 # count node*task configs
    command_list: List[str] = []
    for i in range(len(nodes)):
        n_list: List[int] = nodes[i]
        t_list: List[int] = tasks[i]
        node_str: str = ""
        for n in n_list:
            node_str += str(n) + " "
        tasks_str: str = ""
        for t in t_list:
            tasks_str += str(t) + " "
        command: str = f"./{file} --datapath {datapath} --variant {variant} --size {size} --clusters {clusters} --dimension {dimension} --seed {seed} --nodes {node_str} --ntasks-per-node {tasks_str} --repeat {repeat} --cluster {compute_cluster} --means-set {mean_set} --init-seed {init_seed}"
        print(f"> Running commmand: {command}")
        if not debug:
            subprocess.run(command.split())
        command_list.append(command)
        config_counter += len(n_list) * len(t_list)   
    return config_counter, command_list

def progress(filecount: int, old_results: int) -> None:
    """ Creates a progress bar in terminal to visualize the percentage of jobs started after submitting """
    current: int = 0
    new: int = len(glob.glob1(".","*.out"))
    pbar = tqdm(total=filecount, desc="Waiting for jobs to start")
    time_spend = 0
    while ((new - old_results < filecount) and time_spend < 30*60):
        current = new
        time.sleep(1)
        new = len(glob.glob1(".","*.out"))
        pbar.update(new - current)
        time_spend += 1
    pbar.close()

def sleep_bar(seconds: int, msg="Waiting") -> None:
    """ Sleep for given amount of seconds and visualize with a progress bar """
    for i in tqdm(range(0, seconds), total=seconds, desc=msg):
        time.sleep(1)

def wait_on_queue() -> None:
    """ Wait for all processes of a given account name to exit the queue, or 15 mins.
        Visualized with a progress bar on the 15 min timer. """
    try:
        grep: str = subprocess.check_output(f"squeue | grep {account_name}", shell=True)
    except Exception as e:
        grep: str = ""
    time_spent = 0
    pbar = tqdm(total=15*60, desc="Waiting for jobs to finish")
    while (grep and time_spent < 15*60):
        time.sleep(1)
        time_spent += 1
        pbar.update(1)
        try:
            grep = subprocess.check_output(f"squeue | grep {account_name}", shell=True)
        except Exception as e:
            grep = ""
    pbar.close()
    time.sleep(3)

def write_commands_to_file(command_list: List[str]) -> str:
    """ Writes individual commands (deconstructed) into a file. 
        Used for checking each result individually.
        Can also be used to rerun individual configs manually afterwards. """
    config_list: List[str] = []
    # Construct pairs of command strings and booleans to mark if they have a valid result
    for command in command_list:
        parts: List[str] = command.split("--")
        script: str = parts[0].strip()
        arguments: List[str] = parts[1:]
        names: Dict[str, List[str]] = split_arguments(arguments)
        # Remove arguments with singular value
        template: str = f"{script} --datapath {names['datapath'][0]} --seed {names['seed'][0]} --repeat {names['repeat'][0]} --means-set {names['means-set'][0]} --init-seed {names['init-seed'][0]} --cluster {names['cluster'][0]}"
        del names["datapath"]
        del names["seed"]
        del names["repeat"]
        # Create command for each config
        for variant in names["variant"]:
            for size in names["size"]:
                for clusters in names["clusters"]:
                    for dimension in names["dimension"]:
                        for nodes in names["nodes"]:
                            for tasks in names["ntasks-per-node"]:
                                config: str = f"{template} --variant {variant} --size {size} --clusters {clusters} --dimension {dimension} --nodes {nodes} --ntasks-per-node {tasks}"
                                config_list.append(config)
                            
    # Create file containing all individual commands
    filename = f"commands_{date}.txt"
    with open(filename, "w") as f:
        for c in config_list:
            f.write(c+"\n")
    print(f"> Commands saved in: {filename}")
    return filename

def validate_file(filename: str, repeat=10) -> bool:
    """ Check if the number of retries succeeded in a given *.out filename """
    with open(filename, "r") as f:
        for line in (f.readlines()[-repeat:]):
            if not line.startswith("OK"):
                return False
    return True

def validate_results(filename) -> List[str]:
    """ Searches and validates *.out files corresponding to each command in the given commands.txt file.
        Duplicate and invalid results are removed. """
    invalid: List[str] = []
    if not exists(filename):
        return
    # Check the result of each configuration and report errors
    suffix = "id[0-9]*.out" # output file glob regex
    valid = False
    with open(filename, "r") as f:
        for config in f.readlines():
            a: Dict[str, List[str]] = split_arguments(config.split("--")[1:])
            prefix = f"{a['variant'][0]}_s{a['size'][0]}_c{a['clusters'][0]}_d{a['dimension'][0]}_{a['nodes'][0]}x{a['ntasks-per-node'][0]}_"
            file_list: List[str] = glob.glob1(".",f"{prefix}{suffix}")
            valid = False
            for output_file in file_list:
                if valid:
                    os.remove(output_file) # Other file is already valid, remove duplicate result files
                    print(f"WARNING: duplicate result removed for config:\n {config}")
                elif validate_file(output_file, repeat=int(a['repeat'][0])):
                    valid = True
                else:
                    os.remove(output_file)
                    print(f"WARNING: invalid result removed for config:\n {config}")
            if not valid:
                invalid.append(config.strip("\n"))
    return invalid

def resubmit_jobs(invalid: List[str], filename: str, retries: int = 4) -> List[str]:
    """ Resubmits failed configuration. Limited retries.
        Adds waiting calls, including progress bars. """
    it: int = 0
    while (len(invalid) > 0 and it < retries):
        print(f"> Resubmitting {len(invalid)} jobs")
        for command in invalid:
            subprocess.run(command.split())
        # Wait for jobs to finish up
        old_results = len(glob.glob1(".","*.out"))
        filecount = len(invalid)
        progress(filecount, old_results)
        wait_on_queue()
        invalid = validate_results(filename)
        print(f"{len(invalid)} invalid results encountered")
        it += 1
    return invalid

def sort_results(filename: str) -> None:
    """ Sort results in filename according by execs order.
        Cuts own_inc_loc lines and pastes them at the end
        Only works for original 4 implementations. """
    def convert(s: str):
        try:
            res = int(s)
        except ValueError as e:
            res = float(s)
        return res

    with open(filename, "r") as f:
        lines: List[str] = f.readlines()
        res: List[str] = []
        for e in execs:
            cut: List[Tuple] = [tuple(convert(y) if y not in execs else y for y in x.split()) for x in lines if x.startswith(e+" ")]
            cut = sorted(cut, key=operator.itemgetter(1, 2, 3, 4, 5))
            cut = [" ".join([str(x) for x in t]) for t in cut]
            res = res + ["".join(str(x)) for x in cut]
        content: str = "\n".join(res)

    with open(filename, "w") as f:
        f.write(content)

def process_results(output_dir: str, compute_cluster: str, scriptpath: str, ex_num: int, file_date: str) -> int:
    """ Process the results present in the current directory. 
        Needs the path to process-results.py in scriptpath.
        Results are saved in directory output_dir. """
    if not Path(output_dir).is_dir():
        Path(output_dir).mkdir(parents=True)
    filename = f"{output_dir}/EX{ex_num}-{compute_cluster}-RESULTS-{file_date}.txt"
    with open(filename, "w") as result_file:
        file = scriptpath + "/process-results.py"
        if not exists(file):
            return 1
        subprocess.run([f"./{file}", "."], stdout=result_file)
    print(f"> Running command: ./{file} .")
    sort_results(filename)
    print(f"> Results saved in {filename}")
    return 0

def run_experiment(config: Dict[str, any], options: Dict[str, any], ex_num: int) -> int:
    """ Run experiment ex_num with the provided experiment configuration and command line arguments. """
    print(f"> Running experiment {ex_num} ...")
    
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
    init_seed: int = options["init-seed"]
    file_date: str = options["date"]
    mean_set: int = options["means-set"]

    # Check if any .out files still in directory
    old_results: int = len(glob.glob1(".","*.out"))
    if old_results > 0:
        print(f"WARNING: {old_results} '*.out' files already present in directory!")
    
    # Prepare commands and execute script
    file = scriptpath + "/submit-exp.py"
    if not exists(file):
        return 1
    if len(nodes) != len(tasks) or len(nodes) < 1:
        print(f"ERROR: nodes and tasks should be a list of list of int containing nodes * tasks configurations constructed using the index of the outer list. \
              Example: [[1], [8]] [[2, 4], [3]] refers to 1*2, 1*4, and 8*3", file=sys.stderr)
        return 1
    
    # Create commands to submit jobs for each nodes * tasks config
    config_counter, command_list = submit_jobs(file, datapath, variant, size, clusters, dimension, seed, nodes, tasks, repeat, compute_cluster, mean_set, init_seed)
   
    # Create commands report file
    print("> Writing all executed commands to file ...")
    filename: str = write_commands_to_file(command_list)

    # Wait for jobs to start and finish
    filecount = config_counter * len(variant.split()) * len(size.split()) * len(clusters.split()) * len(dimension.split())
    if not debug:
        progress(filecount, old_results) # waits for all jobs to start
        wait_on_queue() # checks run queue for set account name
    print("> Job runs finished!")
    
    # Validate result of each individual command, rerun failed configs
    print("> Validating command outputs ...")
    invalid: List[str] = validate_results(filename)
    print(f"{len(invalid)} invalid results encountered")
    if not debug:
        invalid = resubmit_jobs(invalid, filename)
        if (len(invalid)):
            s = '\n'.join(invalid)
            print(f"ERROR: not all results are valid. Rerun commands:{s}", file=sys.stderr)
            return 1

    # Process results
    print("> Processing available results ...")
    res: int = process_results(output_dir, compute_cluster, scriptpath, ex_num, file_date)
    if res != 0:
        return res
    
    # Finish up
    print("Done!")
    print(f"Visualize results by running:\n./ex{ex_num}.py --compute-cluster {compute_cluster} --file-date {date} --datapath {output_dir}")
    return 0

def main():
    parser = ArgumentParser()
    parser.add_argument("--debug", "--d", dest="debug", action="store_true",
                            help="Debug mode: no job submits")
    parser.add_argument("--compute-cluster", "--c", dest="compute-cluster", type=str, required=True, choices=["DAS5", "DAS6"],
                            help="Compute cluster in use")
    parser.add_argument("--experiment", "--e", dest="experiment", type=int, default=1, choices=[1, 2],
                            help="Experiment number to execute")
    parser.add_argument("--alternative", "--a", dest="alternative", action="store_true",
                            help="Set alternative configuration for given experiment")
    parser.add_argument("--outputdir", "--o", dest="outputdir", type=str, default="results",
                            help="Location of the resulting output")
    parser.add_argument("--datapath", "--dp", dest="datapath", type=str, default="/var/scratch/dbuurman/kmeans",
                            help="Location of the dataset")
    parser.add_argument("--scriptpath", "--s", dest="scriptpath", type=str, default="../tupl-kmeans-39f1073/support",
                            help="Location to directory of submit-exp.py script")
    parser.add_argument("--clean", dest="clean", action="store_true",
                            help="Clean output scripts; overrides other arguments.")
    parser.add_argument("--process", dest="process", action="store_true",
                            help="Process output scripts; overrides other arguments.")
    parser.add_argument("--validate", dest="validate", action="store_true",
                            help="Only validate results in current folder; overrides other arguments.")
    parser.add_argument("--commands-file", dest="commands-file", type=str, default=f"commands_{date}.txt",
                            help="Clean output scripts; overrides other arguments.")
    parser.add_argument("--init-seed", dest="init-seed", type=int,
                            help="MPI rank init seed")
    parser.add_argument("--means-set", "--m", dest="means-set", type=int, default=1,
                            help="Mean set to use for initialization")
    parser.add_argument("--date", dest="date", type=str, default=date,
                            help="Date used in result file name")
    args = parser.parse_args()
    options: Dict[str, any] = dict(vars(args))

    # Only clean output files
    if options["clean"]:
        for f in glob.glob1(".", "*.out"):
            os.remove(f)
        return 0
    del options["clean"]

    if options["process"]:
        return process_results(options["outputdir"], options["compute-cluster"], options["scriptpath"], options["experiment"], options["date"])
    del options["process"]

    # Only validate results
    if options["validate"]:
        if exists(options["commands-file"]):
            invalid = validate_results(options["commands-file"])
            msg = '\n'.join(invalid)
            print(f"{len(invalid)} invalid results encountered:\n{msg}")
        return 0
    del options["validate"]
    del options["commands-file"]

    # Enable debug mode
    if options["debug"]:
        global debug
        debug = True
        print("--- Debug mode enabled ---")

    # Generate random seed if none set
    if options["init-seed"] == None:
        seed: int = random.randint(0, sys.maxsize)
        print(f"WARNING: --init-seed param not set; random seed {seed} used")
        options["init-seed"] = seed

    # Run selected experiment
    experiment: int = options["experiment"]
    del options["experiment"]
    if experiment == 1:
        if options["alternative"]:
            ex1_config["nodes"] = [[1]] # change to single node comparison experiment
        return(run_experiment(ex1_config, options, 1))
    elif experiment == 2:
        if options["alternative"]:
            ex2_config[options["compute-cluster"]]["size"] = "28" # change to single bigger input size
        return(run_experiment(ex2_config[options["compute-cluster"]], options, 2))

if __name__ == "__main__":
    sys.exit(main())
    