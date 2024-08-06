#!/usr/bin/env python3

# 
# Script to process the TIME EXP output from the Kmeans .out files
# Based on process-results.py by Kristian Rietveld
# 
# Author: Dennis Buurman, Leiden University
import sys
import re
from pathlib import Path

from typing import Dict

# Taken from common.py by Kristian Rietveld
def parse_job_name(jobname : str) -> dict:
    m = re.match(r"([\w_]+)_s(\d+)_c(\d+)_d(\d+)_(\d+)x(\d+)_id.*", jobname)
    if not m:
        raise ValueError("Cannot parse jobname '{}'.".format(jobname))

    labels = ["variant", "size", "clusters", "dimension", "n_nodes",
              "ntasks_per_node"]
    params = { k: v for k, v in zip(labels, m.groups()) }
    return params

def process_file(filename: Path) -> None:
    print("Processing file:", filename, file=sys.stderr)

    params: Dict = parse_job_name(filename.name)
    
    results: Dict[str, any] = dict()
    # Read lines
    with filename.open() as f:
        for line in f:
            # Extract values
            if line.startswith("EXP "):
                line = line.rstrip()

                m = re.match(r"EXP (\d+): ", line)
                if not m:
                    continue

                exp = int(m.group(1))
                if not exp in results:
                    results[exp] = dict()

                # Strip match EXP (\d+) part
                line = line[m.span(0)[1] : ]

                m = re.fullmatch(r"iterations (\d+)", line)
                if m:
                    results[exp]["iterations"] = int(m.group(1))

            elif line.startswith("TIME EXP"):
                line = line.rstrip()
                m = re.match(r"TIME EXP (\d+): ", line)
                if not m:
                    continue
                
                exp = int(m.group(1))
                if exp not in results:
                    results[exp] = dict()

                line = line.split(':')[1].strip()
                m = re.fullmatch(r"\[.+\] (\w+) (\d+\.\d+) s", line)
                if m:
                    if m.group(1) not in results[exp]:
                        results[exp][m.group(1)] = list()
                    results[exp][m.group(1)].append(float(m.group(2)))
                
                m = re.fullmatch(r"\[.+\] (\w+) (\d+\.(\d+|\d+e-\d+)) s \| (\d+\.(\d+|\d+e-\d+)) s averaged", line)
                if m:
                    if m.group(1) not in results[exp]:
                        results[exp][m.group(1)] = list()
                    results[exp][m.group(1)].append(float(m.group(2)))
                    if str(m.group(1))+"Avg" not in results[exp]:
                        results[exp][str(m.group(1))+"Avg"] = list()
                    results[exp][str(m.group(1))+"Avg"].append(float(m.group(4)))

    # Write output
    # variant, size, clusters, dimension, n_nodes, ntasks_per_node, 
    # init, calc, reassign, communication, recalc, iterations

    for exp in results:
        worstcalctime = max(results[exp]["calculationTime"])
        # worstcalctimeavg = max(results[exp]["calculationTimeAvg"])
        worstreassigntime = max(results[exp]["reassignTime"])
        # worstreassigntimeavg = max(results[exp]["reassignTimeAvg"])
        worstcommunicationtime = max(results[exp]["communicationTime"])
        # worstcommunicationtimeavg = max(results[exp]["communicationTimeAvg"])
        worstmeanrecalctime = max(results[exp]["meanRecalculationTime"])
        # worstmeanrecalctimeavg = max(results[exp]["meanRecalculationTimeAvg"])

        print(params["variant"], params["size"], params["clusters"],
              params["dimension"], params["n_nodes"], params["ntasks_per_node"],
              max(results[exp]["initTime"]), worstcalctime, worstreassigntime, worstcommunicationtime,
              worstmeanrecalctime, results[exp]["iterations"])

def main() -> int:
    if len(sys.argv) != 2:
        print("usage: {} <results dir>".format(sys.argv[0]), file=sys.stderr)
        return 1
    
    results_dir: Path = Path(sys.argv[1])

    print("variant, size, clusters, dimension, n_nodes, ntasks_per_node, init_time, calculation_time, reassign_time, communication_time, recalculation_time, iterations")
    for filename in results_dir.glob("*.out"):
        process_file(filename)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())