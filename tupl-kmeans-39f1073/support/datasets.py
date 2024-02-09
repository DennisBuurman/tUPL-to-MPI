#!/usr/bin/env python3

#
# Author: Kristian Rietveld, Leiden University
#

import os
from pathlib import Path
import common
from common import DataSetRegistry

import sys
from argparse import ArgumentParser

# generator_exec = Path("..") / "generator" / "generator"
generator_exec = Path("../../dataset_generator/generator")


def print_query(params):
    print("Query: " + DataSetRegistry.summarize_query(params))

def print_matches(datapath : Path, matches):
    print("The following datasets match:")
    for path, properties in matches:
        path = path.relative_to(datapath)
        print("    {} (.../{})".format(DataSetRegistry.summarize_query(properties), path))

def list_datasets(datapath : Path, params : dict) -> int:
    registry = DataSetRegistry(datapath)

    print_query(params)
    print_matches(datapath, registry.query(params))

    return 0


def create_datasets(datapath : Path, params : dict) -> int:
    # All required options must be present in params
    if not all(p in params for p in common.DataSetRegistry.required_properties):
        print("Error: All properties must be specified to generate datasets.", file=sys.stderr)
        return 1

    if not generator_exec.exists():
        print("Error: generator executable '{}' does not exist.".format(generator_exec), file=sys.stderr)
        return 1

    registry = DataSetRegistry(datapath)

    # Iterate combinations of parameters and generate the necessary
    # datasets
    for p in common.iter_param_combinations(params):
        print("Considering " + common.DataSetRegistry.summarize_query(p))
        if registry.find_datafiles(p):
            print("  Data file already present, skipping generation.")
        else:
            print("  Generating data file ...")

            dirname = datapath / DataSetRegistry.get_dir_name(p)
            print("    + mkdir -p {}".format(dirname))
            dirname.mkdir(parents=True)

            cmd = "{exec} {seed} {size} {clusters} {dimension} {outdir}".format(exec=generator_exec.absolute().as_posix(), seed=p["seed"], size=p["size"], clusters=p["clusters"], dimension=p["dimension"], outdir=dirname)
            print("    + " + cmd)
            os.system(cmd)

    return 0


def main() -> int:
    # Set up argument parser
    parser = ArgumentParser()
    DataSetRegistry.add_dataset_parameters(parser)

    parser.add_argument('command', metavar='command', type=str,
                        choices=['list', 'create'],
                        help='command to execute (list or create)')

    args = parser.parse_args()

    # Check mandatory datapath argument
    if not args.datapath:
        print("Error: no datapath specified.", file=sys.stderr)
        return 1

    # Create params dictionary that only contains dataset parameters
    datapath = Path(args.datapath)
    command = str(args.command)

    params = { key: val for key, val in vars(args).items() if val != None }
    del params["datapath"]
    del params["command"]

    # Run desired command
    if command == "create":
        return create_datasets(datapath, params)
    elif command == "list":
        return list_datasets(datapath, params)

    return 1


if __name__ == '__main__':
    sys.exit(main())
