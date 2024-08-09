This folder contains the code to run the benchmarks.
These benchmarks contain the first two experiments from [1] and new experiments.

Quick start:
# Change account_name variable in benchmark.py to your account_name
mkdir results
chmod +x benchmark.py ex1.py ex2.py
./benchmark.py -c DAS5 -e 1
chmod +x input-size-variation.py thread-count-variation.py timing-stack.py

Default parameters of benchmark.py:
# 5 iteration run:
./benchmark.py --c <cluster> --e <ex_num> --init-seed 340632450 --m 3
# ~70 iteration run:
./benchmark.py --c <cluster> --e <ex_num> --init-seed 340632450 --m 7
# <cluster> can be either DAS5 or DAS6
# <ex_num> options can be found in benchmark.py

Visualization:
./input-size-variation.py   --compute-cluster <cluster> --file-date <date> --datapath <path> --ex-num <ex_num>
./thread-count-variation.py --compute-cluster <cluster> --file-date <date> --datapath <path> --ex-num <ex_num>
./timing-stack.py           --compute-cluster <cluster> --file-date <date> --datapath <path> --ex-num <ex_num>
# <cluster> can be either DAS5 or DAS6
# <date> date suffix from results file produced by benchmark.py (format 'dd-mm-yyyy')
# <ex_num> options can be found in benchmark.py

After running benchmark.py, commands for visualizing the produced results are outputted. 
Execute the printed command(s) to create the graphs corresponding to the experiment performed.

The benchmark.py script uses the submit-exp.py script from the tupl-kmeans directory (created by A. Hommelberg, K.F.D. Rietveld).
Adding implementation variants to benchmark.py therefore requires you to also add them to submit-exp.py!

[1] A. Hommelberg, B. van Strien, K. F. D. Rietveld, and H. A. G. Wijshoff, “A new framework for expressing, parallelizing and optimizing big data applications,” arXiv preprint arXiv:2203.01081, 2022.