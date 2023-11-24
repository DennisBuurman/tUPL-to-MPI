This is a fork of the original repository by Anne Hommelberg in order
to clean things up. The code clean up was done by Kristian Rietveld.
Also, several support scripts have been written in Python to aid in
generating data sets, submitting jobs and processing job logs.

These programs are meant to be executed on MPI-based clusters and have
been tested on DAS-4.


## How to run an experiment

Start by compiling the programs in the `generator` and `src` directories.
On DAS you need to have an "openmpi" module loaded and perhaps a gcc more
recent than the system gcc as well:

```
module load gcc/6.4.0
module load openmpi/gcc
```

### Data set generation

Next, we need to generate data files. On DAS, you want to store this in
`/var/scratch/accountname/kmeans`. To generate data a useful script is
included in the `support` directory. Let's generate a number of data files
for size 2^24 and different amounts of clusters:

```
mkdir /var/scratch/accoutname/kmeans
cd support
./datasets.py create --datapath /var/scratch/accountname/kmeans \
    --size 24 --clusters 4 8 12 24 32 --dimension 4 --seed 771
```

This will take a while and results in ~4 GiB of data. (Note that larger
sizes will result in much larger data sets!). Let's list the currently
available data sets:

```
./datasets.py list --datapath /var/scratch/accountname/kmeans
```

We can even add parameters to query a subset of the available datasets:

```
./datasets.py list --datapath /var/scratch/accountname/kmeans --clusters 4 8
```

### Running experiments

The `submit-exp.py` script can generate job scripts and automatically submit
these to the job scheduler. When multiple experiment parameters are
specified, multiple jobs will be submitted. For instance, for the data sets
that we have just created:

```
./submit-exp.py --datapath /var/scratch/accountname/kmeans \
    --variant own --repeat 5 --size 24 --clusters 4 8 12 24 32 --seed 771
```

The script will verify whether the required data sets exist and if not it
will tell you. The jobs will be submitted immediately. If you want to have
the script write the jobs to files instead, add the `-o` option to specify
an output directory for the job scripts. The job scripts can then be
inspected with an editor and submitted at a later time using the `sbatch`
command on DAS.

Files with the suffix ".out" are created by the jobs. You likely want to
put files belonging to a single experiment in a directory. Say `results`.
A simple script is supplied to collect and summarize results from a series
of job logs in a directory:

```
./process-results.py results
```

This summary can then be imported into a spreadsheet or otherwise to be
aggregated and to create plots. Further modifications to
`process-results.py` might be necessary to improve the output.
