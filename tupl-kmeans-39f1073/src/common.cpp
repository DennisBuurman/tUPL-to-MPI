/*
 * MPI K-means implementation
 *
 * Author: Anne Hommelberg
 *
 * Code refactoring done by Kristian Rietveld, Leiden University.
 */

#include "common.h"

#include "mpi-utils.h"

#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <getopt.h>


int mpi_rank, mpi_size;
char mpi_hostname[MPI_MAX_PROCESSOR_NAME];

MPI_Datatype MPI_DATAPOINT;


double ** allocate2dDoubleArray(uint64_t rows, int columns) {
  double * data = new double[rows * columns];
  double ** array = new double*[rows];
  for (uint64_t i = 0; i < rows; i++) {
    array[i] = &(data[columns*i]);
  }
  return array;
}

void initRandom(const struct Options &options)
{
  unsigned int seed = 0;
  if (mpi_rank == 0) {
    if (!options.seedFlag) {
      std::ifstream randomIn("/dev/random");
      randomIn.read(reinterpret_cast<char *>(&seed), sizeof(seed));
      randomIn.close();
    } else {
      seed = options.seed;
    }

    std::cout << "EXP " << options.currentRun << ": using seed "
        << std::hex << seed << std::dec << std::endl;
  }

  MPI_Bcast(&seed, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

  /* All MPI processes initialize using the same seed */
  srand(seed);
}

static void
showHelp(const char *progName)
{
  if (mpi_rank == 0)
    std::cerr << "Usage: " << progName << " "
        << "-i input_dir -f format -k number_of_means -d convergence_delta "
        << "-t threshold -s exp_suffix [-r #runs]"
        << std::endl;
}

/* We need this enum to have an easy way to determine whether all
 * options are supplied (all are required).
 */
enum Args
{
  ARGS_NONE = 0,
  ARGS_INPUT_DIR = 1,
  ARGS_CHINESE_FORMAT = 1 << 1,
  ARGS_NUM_MEANS = 1 << 2,
  ARGS_CONVERGENCE_DELTA = 1 << 3,
  ARGS_THRESHOLD_MULTIPLIER = 1 << 4,
  ARGS_OUTPUT_SUFFIX = 1 << 5,

  ARGS_ALL = (ARGS_OUTPUT_SUFFIX << 1) - 1
};

bool parseArgs(int argc, char ** argv, struct Options &options)
{
  char c;
  enum Args providedArgs = ARGS_NONE;
  const char *progName = argv[0];

  options.numRuns = 1;

  while ((c = getopt(argc, argv, "hi:f:k:d:t:s:r:m:x:")) != -1)
    {
      switch (c)
        {
          case 'i':
            options.inputDir = optarg;
            providedArgs = (Args)(providedArgs | ARGS_INPUT_DIR);
            break;

          case 'f':
            options.chineseFormat = (std::atoi(optarg) != 0);
            providedArgs = (Args)(providedArgs | ARGS_CHINESE_FORMAT);
            break;

          case 'k':
            options.numMeans = std::atoi(optarg);
            providedArgs = (Args)(providedArgs | ARGS_NUM_MEANS);
            break;

          case 'd':
            options.convergenceDelta = std::atof(optarg);
            providedArgs = (Args)(providedArgs | ARGS_CONVERGENCE_DELTA);
            break;

          case 't':
            options.thresholdMultiplier = std::atof(optarg);
            providedArgs = (Args)(providedArgs | ARGS_THRESHOLD_MULTIPLIER);
            break;

          case 's':
            options.outputSuffix = optarg;
            providedArgs = (Args)(providedArgs | ARGS_OUTPUT_SUFFIX);
            break;

          case 'r':
            options.numRuns = std::atoi(optarg);
            break;

          case 'x':
            options.seedFlag = true;
            options.seed = std::atoi(optarg);
            break;

          case 'm':
            options.meansFlag = true;
            options.meansSet = std::atoi(optarg);
            break;

          case 'h':
            showHelp(progName);
            exit(0);
        }
    }
  
  argc -= optind;
  argv += optind;

  if (argc > 1) {
    if (mpi_rank == 0)
      std::cerr << "Error: not all arguments were parsed." << std::endl;

    return false;
  }

  if (providedArgs != ARGS_ALL) {
    if (mpi_rank == 0)
      std::cerr << "Error: some arguments are missing, experiment configuration not complete." << std::endl;
    showHelp(progName);

    return false;
  }

  if (mpi_rank == 0) {
    std::cout << "DESC: dir=" << options.inputDir << " "
              << "numMeans=" << options.numMeans << " "
              << "convergenceDelta=" << options.convergenceDelta << " "
              << "thresholdMultiplier=" << options.thresholdMultiplier
              << std::endl;
    
    std::cout << "INITIALIZATION: ";
    if (options.meansFlag) {
      std::cout << "MEAN SET " << options.meansSet;
    } else {
      std::cout << "RANDOM";
    }
    std::cout << std::endl;
  }

  return true;
}

int runVariant(int argc, char **argv, kMeansVariantFunc variantFunc)
{
  MPI_Init(&argc, &argv);
  
  int hostNameLength;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Get_processor_name(mpi_hostname, &hostNameLength);
  mpi_hostname[hostNameLength] = 0;
  
  struct Options options;
  // options.seedFlag = false;
  // options.seed = 0;
  options.meansFlag = false;
  options.meansSet = 0;
  if (!parseArgs(argc, argv, options))
    return 1;
  
  double start_time, end_time; 
  
  //some initial output
  start_time = MPI_Wtime();
  std::cout << printPreamble() << " Running on " << mpi_hostname << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) 
    std::cout << std::endl;
  
  //run the algorithm
  int rc = variantFunc(options);
  
  end_time = MPI_Wtime();
  std::cout << "EXP: " << printPreamble() << " fullTime "
            << end_time - start_time << " seconds" << std::endl;
  
  MPI_Finalize();
  return rc;
}


void initFile(struct Options &options, struct FileState &state)
{
  std::stringstream inputFile;
  inputFile << options.inputDir << "/data.txt";
  std::cout << "Reading file: " << inputFile.str() << std::endl;

  options.numDataPoints = 0;

  state.start_time = MPI_Wtime();
  state.in.open(inputFile.str());
  if (!state.in.is_open()) 
    throw std::runtime_error("Could not open input file.");

  //use the first line to determine the dimension of the data;
  std::string line;
  std::getline(state.in,line);
  options.dataDim = std::count(line.begin(), line.end(), '.');
  state.in.seekg(0);

  if (options.dataDim < DATA_DIM) {
    std::cerr << "Warning: data dimension lower than expected, not utulizing all allocated storage space." << std::endl;
  }
  if (options.dataDim > DATA_DIM) {
    throw std::runtime_error("Allocated storage space insufficient, data dimension too high");
  }
}

void finishFile(struct Options &options, struct FileState &state)
{
  double end_time = MPI_Wtime();
  std::cout << "EXP: readTime " << end_time - state.start_time << " seconds" << std::endl;
  std::cout << "Detected " << options.numDataPoints << " data points of dimension " << options.dataDim << "." << std::endl << std::endl;
}

static void bcastDataSize(struct Options &options)
{
  MPI_Bcast(&options.numDataPoints, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&options.dataDim, 1, MPI_INT, 0, MPI_COMM_WORLD);

  uint64_t from, to;
  std::tie(from,to) = mpi_range(options.numDataPoints);
  options.numLocalDataPoints = to - from;
}

static void scatterData(struct Options &options,
                        const void *sendbuf, void *recvbuf, MPI_Datatype mpitype,
                        const size_t sendCountInterval,
                        const size_t displacementInterval)
{
  uint64_t from, to;
  std::tie(from,to) = mpi_range(options.numDataPoints);

  //communicate data to all processes
  int sendCounts[mpi_size];
  int displacements[mpi_size];
  for (int p = 0; p < mpi_size; p++) {
    uint64_t fromp, top;
    std::tie(fromp, top) = mpi_range(options.numDataPoints, p);
    sendCounts[p] = (top - fromp) * sendCountInterval;
    displacements[p] = from * displacementInterval;
  }
  //NOTE: should be adapted to handle more than INT_MAX at a time
  MPI_Scatterv(sendbuf, sendCounts, displacements, mpitype, recvbuf,
               options.numLocalDataPoints * options.dataDim, mpitype,
               0, MPI_COMM_WORLD);
}


double ** readDataAsDoubles(struct Options &options) {
  MPI_Barrier(MPI_COMM_WORLD); //NOTE necessary for alligning output
  
  //have the first process read all data
  std::vector<double> alldata;
  if (mpi_rank == 0) {
    struct FileState state;
    initFile(options, state);

    double buf[options.dataDim];
    while (readPoint(options, state, &buf[0])) {
      for (int i = 0; i < options.dataDim; ++i)
        alldata.emplace_back(buf[i]);
    }

    finishFile(options, state);
  }
  
  bcastDataSize(options);

  //allocate the local data
  double ** data = allocate2dDoubleArray(options.numLocalDataPoints, options.dataDim);

  scatterData(options, alldata.data(), data[0], MPI_DOUBLE,
              options.dataDim, options.dataDim);
  
  if (mpi_rank == 0) {
    alldata.clear();
    std::cout << std::endl;
  }
  return data;
}


void defineMPIDataPoint(void) {
  //prepare a MPI_Datatype for the DataPoint struct
  MPI_Datatype type[2] = {MPI_DOUBLE, MPI_INT};
  int blockLength[2] = {DATA_DIM, 1};
  MPI_Aint displacement[2] = {0, DATA_DIM*sizeof(double)};
  MPI_Type_create_struct(2, blockLength, displacement, type, &MPI_DATAPOINT);
  MPI_Type_commit(&MPI_DATAPOINT);
}

DataPoint * readDataAsDataPoints(struct Options &options) {
  MPI_Barrier(MPI_COMM_WORLD); //NOTE necessary for alligning output
  
  //have the first process read all data
  std::vector<DataPoint> alldata;
  if (mpi_rank == 0) {
    struct FileState state;
    initFile(options, state);

    DataPoint element;
    while (readPoint(options, state, element.vector))
      alldata.push_back(element);

    finishFile(options, state);
  }
  
  bcastDataSize(options);

  //allocate and receive the local data
  DataPoint * data = new DataPoint[options.numLocalDataPoints];

  scatterData(options, alldata.data(), data, MPI_DATAPOINT,
              1, sizeof(DataPoint));
  
  if (mpi_rank == 0) {
    alldata.clear();
    std::cout << std::endl;
  }
  return data;
}



//returns the squared Eucledian distance, note that to compare distances no square root is necessary
double calculateDistance(const double * vector1, const double * vector2,
                         const int size) {
  double distance = 0;
  for (int i = 0; i < size; i++) {
    distance += ((vector1[i] - vector2[i]) * (vector1[i] - vector2[i]));
  }
  return distance;
}
