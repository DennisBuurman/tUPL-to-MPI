#!/usr/bin/env python3

#
# Author: Kristian Rietveld, Leiden University
#

from pathlib import Path
import sys
from string import Template
import subprocess
from argparse import ArgumentParser

import common

defaults = { "size": 28, "clusters": 4, "dimension": 4, "seed": 971 }

modules = {
"DAS5" : """. /etc/bashrc
. /etc/profile.d/modules.sh
module load openmpi/gcc/64
""",
"DAS6" : """. /etc/bashrc
. /etc/profile.d/lmod.sh
module load openmpi/gcc/64
MPI_RUN=mpirun
"""
}

job_template = """#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH -N $n_nodes
#SBATCH --ntasks-per-node=$ntasks_per_node
#SBATCH -J kmeans_$jobname
#SBATCH -o ${jobname}_id%j.out

$config

$modules

datadir=$datadir
numRuns=$repeat

APP=${script_path}/../src/$appexec
## ./MPI_Kmeans -i input_file -f format -k number_of_means -d convergence_delta -t threshold -s exp_suffix -r numRuns -x init_seed
ARGS="-i $$datadir -f 1 -k $clusters -d 0.0001 -t 0 -s ${suffix} -r $$numRuns -m $means_set -x $init_seed"
OMPI_OPTS="--mca btl ^usnic,tcp"

$$MPI_RUN $$OMPI_OPTS $$APP $$ARGS

if [ $$? != 0 ]; then
    echo "ERR: job failed";
    exit 1
fi

echo
${script_path}/validate-output.py $$datadir "${variant}${suffix}" $$numRuns
"""


execs = {
    "own": "MPI_Kmeans",
    "own_loc": "MPI_Kmeans_localized",
    "own_inc": "MPI_Kmeans_incremental",
    "own_inc_loc": "MPI_Kmeans_incremental_localized",
    "own_values_only": "MPI_Kmeans_values_only",
    "own_im": "MPI_Kmeans_no_updates",
    "own_m": "MPI_Kmeans_local_values"
}



def submit_job(job_script : str):
    '''Submit given job_script as job.'''
    process = subprocess.Popen("sbatch", stdin=subprocess.PIPE)
    process.stdin.write(bytes(job_script, 'UTF-8'))
    process.stdin.close()
    if process.wait() != 0:
        print("error submitting job.")
        exit(1)

def dump_config(params : dict) -> str:
    s = "# Experiment configuration:\n"
    for key, value in params.items():
        s += "#    {}={}\n".format(key, value)
    return s

def create_job_script(params) -> str:
    '''Create a job script using the given parameters.'''
    global job_template, execs

    templ = Template(job_template)
    job_script = templ.substitute(config=dump_config(params),
                                  modules=modules[params["cluster"]],
                                  jobname=common.create_job_name(params),
                                  script_path=Path(__file__).absolute().parent.as_posix(),
                                  appexec=execs[params["variant"]],
                                  suffix="test",   # TODO: make option?
                                  init_seed=params["init-seed"],
                                  means_set=params["means-set"],
                                  **params)
    return job_script


def enumerate_job_scripts(datapath : Path, options : dict) -> list:
    '''Given a dictionary of experiment options, enumerate
    all possible combinations of options and generate job scripts
    for these.'''

    # Apply defaults if option was not set
    for key, val in defaults.items():
        if key not in options or options[key] == None:
            options[key] = val

    # Ensure every option value is a list (even if it contains a single value)
    for key, val in options.items():
        if not isinstance(val, list):
            options[key] = [val]

    # Now generate jobscripts using datafiles form the registry
    all_ok = True
    registry = common.DataSetRegistry(datapath)
    job_scripts = list()

    # Enumerate every possible combination of the options to generate
    # the experiments to be performed.
    for params in common.iter_param_combinations(options):
        datafiles = registry.find_datafiles(params)
        if not datafiles:
            print("Could not find datafile for configuration: " + common.DataSetRegistry.summarize_query(params))
            all_ok = False
        else:
            for datadir in datafiles:
                if not datadir.exists():
                    print("Data directory does not exist: " + datadir)
                    all_ok = False
                if not (datadir / "data.txt").exists():
                    print("Data file in data directory does not exist (?): " + datadir)
                    all_ok = False

                params["datadir"] = datadir.absolute().as_posix()
                job_scripts.append((common.create_job_name(params),
                                   create_job_script(params)))

    # Only return success (and have caller proceed) if all job scripts
    # were generated successfully.
    if all_ok:
        return job_scripts
    return []


def generate_job_scripts(datapath : Path, output_dir : str, options : dict) -> int:
    if output_dir:
        if not Path(output_dir).is_dir():
            print("Error: output directory '{}' is not a directory.".format(output_dir), file=sys.stderr)
            return 1

    job_scripts = enumerate_job_scripts(datapath, options)

    if output_dir:
        print("Writing job scripts to directory '{}' ...".format(output_dir))

        for name, script in job_scripts:
            filename = Path(output_dir) / (name + ".job")
            print(filename)
            with filename.open("w") as fh:
                fh.write(script)
    else:
        print("Submitting jobs to scheduler ...".format(output_dir))

        for name, script in job_scripts:
            submit_job(script)

    return 0


def main() -> int:
    parser = ArgumentParser()
    parser.add_argument("-o", dest="output_dir",
                        help="Output directory for job scripts, if not specified jobs will be submitted immediately.")

    parser.add_argument("--nodes", dest="n_nodes",
                        nargs="*", default="8",
                        help="Number of nodes")
    parser.add_argument("--ntasks-per-node", dest="ntasks_per_node",
                        nargs="*", default="8",
                        help="Number tasks to start on each node")
    parser.add_argument("--variant", dest="variant",
                        nargs="*", default="own", choices=execs,
                        help="Variant of the algorithm to run")
    parser.add_argument("--repeat", dest="repeat", type=int,
                        nargs=1, default=5,
                        help="Number of times to repeat the experiment")
    parser.add_argument("--cluster", dest="cluster", type=str,
                        nargs=1, default="DAS5", choices=["DAS5", "DAS6"],
                        help="Cluster the experiment is run on")
    parser.add_argument("--init-seed", dest="init-seed", type=int,
                        nargs=1, default=340632450,
                        help="Initialization seed for the MPI ranks")
    parser.add_argument("--means-set", dest="means-set", type=int, nargs=1,
                        help="Mean set to use as initialization")
    common.DataSetRegistry.add_dataset_parameters(parser)

    # Parse arguments & enumerate job script
    args = parser.parse_args()
    if not args.datapath:
        print("Error: no datapath specified.", file=sys.stderr)
        return 1

    # Prepare options for the generation of job scripts by removing
    # arguments that are no experiment parameters.
    options = dict(vars(args))

    # Remove datapath and output_dir options because these are separate from
    # all other options and do not need to be enumerated.
    datapath = Path(options['datapath'])
    del options['datapath']

    output_dir = options['output_dir']
    del options['output_dir']


    return generate_job_scripts(datapath, output_dir, options)

if __name__ == "__main__":
    sys.exit(main())
