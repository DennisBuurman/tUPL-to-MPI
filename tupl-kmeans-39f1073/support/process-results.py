#!/usr/bin/env python3

#
# A simple script to extract/summarize results from a directory containing
# job logs. The output is similar to Anne's process_results.cpp, however
# more data columns have been added and the way "fullTime" is determined
# differs.
#
#
# Author: Kristian Rietveld, Leiden University
#

import sys
import re
from pathlib import Path

import common

def process_file(filename : Path) -> None:
    print("Processing file:", filename, file=sys.stderr)

    params = common.parse_job_name(filename.name)

    global_results = dict()
    exp_results = dict()    # type: dict

    # Process all lines
    with filename.open() as fh:
        for line in fh:
            if line.startswith("EXP:"):
                line = line.rstrip()
                m = re.fullmatch(r"EXP: (\w+) (\d+\.\d+) seconds", line)
                if m:
                    global_results[m.group(1)] = float(m.group(2))

                m = re.fullmatch(r"EXP: \[.+\] (\w+) (\d+\.\d+) seconds", line)
                if m:
                    if m.group(1) not in global_results:
                        global_results[m.group(1)] = list()

                    global_results[m.group(1)].append(float(m.group(2)))

            elif line.startswith("EXP "):
                line = line.rstrip()

                m = re.match(r"EXP (\d+): ", line)
                if not m:
                    continue

                exp = int(m.group(1))
                if not exp in exp_results:
                    exp_results[exp] = dict()

                # Strip match EXP (\d+) part
                line = line[m.span(0)[1] : ]

                m = re.fullmatch(r"using seed ([\d\w]+)", line)
                if m:
                    exp_results[exp]["seed"] = m.group(1)
                m = re.fullmatch(r"iterations (\d+)", line)
                if m:
                    exp_results[exp]["iterations"] = int(m.group(1))
                m = re.fullmatch(r"writeTime (\d+.\d+) seconds", line)
                if m:
                    exp_results[exp]["writeTime"] = float(m.group(1))
                m = re.fullmatch(r"\[.+\] (\w+) (\d+\.\d+) seconds", line)
                if m:
                    if m.group(1) not in exp_results[exp]:
                        exp_results[exp][m.group(1)] = list()

                    exp_results[exp][m.group(1)].append(float(m.group(2)))
            elif line.startswith("ERR "):
                line = line.rstrip()
                m = re.fullmatch(r"ERR (\d+): result does NOT match reference", line)
                if m:
                    exp_results[int(m.group(1))]["notMatch"] = True

    # Write output line for each experiment
    for exp in exp_results.keys():
        error = ""
        if len(exp_results[exp]["calculationTime"]) != len(global_results["fullTime"]):
            print("Mismatch number of calculationTime and fullTime reported for experiment {}, marking as error.".format(exp), file=sys.stderr)
            error = "error"

        if "notMatch" in exp_results[exp]:
            error += " incorrect-answer"

        worstcalctime = max(exp_results[exp]["calculationTime"])

        # We define the execution time as read time + calculation time
        # + write time. This is a bit different from how Anne handled
        # things. This is caused by the fact that we allow repeated runs
        # within the MPI executable (so re-using data read into memory).

        worstexectime = global_results["readTime"] + worstcalctime + exp_results[exp]["writeTime"]

        # TODO: FIXME: perhaps we also want to report avg and stdev,
        # instead of just the "worst time"?
        print(params["variant"], params["size"], params["clusters"],
              params["dimension"], params["n_nodes"], params["ntasks_per_node"],
              worstcalctime, worstexectime,
              exp_results[exp]["iterations"],
              error)


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: {} <results dir>".format(sys.argv[0]), file=sys.stderr)
        return 1

    results_dir = Path(sys.argv[1])
    if not results_dir.exists() or not results_dir.is_dir():
        print("error: {} does not exist or is not a directory.".format(results_dir), file=sys.stderr)
        return 1

    for filename in results_dir.glob("*.out"):
        process_file(filename)

    return 0


if __name__ == '__main__':
    sys.exit(main())
