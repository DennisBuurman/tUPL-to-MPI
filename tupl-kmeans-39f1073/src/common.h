/*
 * MPI K-means implementation
 *
 * Author: Anne Hommelberg
 *
 * Code refactoring done by Kristian Rietveld, Leiden University.
 */

#ifndef __COMMON_H__
#define __COMMON_H__

#include <vector>
#include <string>
#include <functional>
#include <fstream>
#include <sstream>

#include <cstdlib>
#include <cstdint>
#include <cassert>

#include <iostream>
#include <cstring>

#include "mpi-utils.h"


/* Common helper functions */

double **allocate2dDoubleArray(uint64_t rows, int columns);
double calculateDistance(const double * vector1, const double * vector2,
                         const int size);
void initRandom(const struct Options &options);

/* Parse command line arguments */

struct Options {
  const char *inputDir;
  bool chineseFormat;
  int numMeans;
  double convergenceDelta;
  double thresholdMultiplier;
  const char *outputSuffix;

  int numRuns;
  int currentRun;

  bool seedFlag;
  int seed;
  bool meansFlag;
  int meansSet;

  /* Determined when reading input file */
  uint64_t numDataPoints;
  uint64_t numLocalDataPoints;
  int dataDim;
};

bool parseArgs(int argc, char ** argv, struct Options &options);

using kMeansVariantFunc = std::function<int(struct Options &options)>;
int runVariant(int argc, char **argv, kMeansVariantFunc variantFunc);


/*
 * Support for non-localized data structure
 */

/* Read data as sequence of doubles (non-localized) */
double ** readDataAsDoubles(struct Options &options);

/* Accessor functions to be inlined in the common templates below. */

static inline double getDataPoint(double **data, const int i, const int d)
{
  return data[i][d];
}

static inline double * getData(double **data, const int i)
{
  return data[i];
}

/* FIXME: shouldn't i be of type uint64_t? */
static inline void setMean(int *belongsToMean, const int i, const int mean)
{
  belongsToMean[i] = mean;
}

static inline int getMean(const int *belongsToMean, const int i)
{
  return belongsToMean[i];
}

/*
 * Support for localized data structure
 */

const int DATA_DIM = 4;

struct DataPoint {
  double vector[DATA_DIM];
  int mean;
};

void defineMPIDataPoint(void);

/* Read data as sequence of DataPoints (localized) */
DataPoint * readDataAsDataPoints(struct Options &options);

/* Accessor functions to be inlined in the common templates below. */

static inline double getDataPoint(DataPoint *data, const int i, const int d)
{
  return data[i].vector[d];
}

static inline double * getData(DataPoint *data, const int i)
{
  return data[i].vector;
}

static inline void setMean(DataPoint *data, const int i, const int mean)
{
  data[i].mean = mean;
}

static inline int getMean(const DataPoint *data, const int i)
{
  return data[i].mean;
}


/*
 * Common I/O functions
 */

struct FileState {
  double start_time;
  std::ifstream in;
};

void initFile(struct Options &options, struct FileState &state);
void finishFile(struct Options &options, struct FileState &state);

template <typename T>
bool readPoint(struct Options &options, struct FileState &state, T dst) {
  if (options.chineseFormat) {
    const int bufSize(128);
    char buf[bufSize];

    if (!state.in.getline(buf, bufSize))
      return false;

    char *start = buf, *pos = buf;

    /* Ignore first element, which is the data point sequence number */
    while (*pos != ' ')
      ++pos;
    *pos = 0;

    ++pos;
    start = pos;

    /* Read point coordinates, one by one. */
    int dim = 0;
    while (*pos != 0) {
      if (*pos == ' ' || *pos == '\n') {
        *pos = 0;
        dst[dim] = std::atof(start);
        ++dim;

        ++pos;
        start = pos;
      } else
        ++pos;
    }

    assert(dim == options.dataDim);
    ++options.numDataPoints;
  }
  else {
    double element;

    if (!(state.in >> element))
      return false;
    options.numDataPoints++;
    dst[0] = element;
      
    for (int d = 1; d < options.dataDim; d++) {
      state.in >> element;
      dst[d] = element;
    }
  }

  return true;
}

template <typename B>
void writeData(const struct Options &options, const std::string variant,
               double ** meanValues,
               const uint64_t * meanSize,
               const B * belongsToMean,
               const bool writeMembership = false)
{
  double start_time = MPI_Wtime();

  //print output to file
  //TODO use MPI_File_write_ordered instead? (setting up a big char array with ascii message may take too much memory to write in one go?)
  if (mpi_rank == 0) {
    std::stringstream ss;
    ss << options.inputDir << "/cluster_centres_" << variant << options.outputSuffix;
    if (options.numRuns != 1)
      ss << "_" << options.currentRun;
    std::ofstream out(ss.str());
    if (!out.is_open()) 
      throw std::runtime_error("Could not open output file for cluster centers.");
    
    for (int i = 0; i < options.numMeans; i++) {
      out << i << " ";
      for (int d = 0; d < options.dataDim; d++) {
        out << meanValues[i][d] << " ";
      }
      out << ", size = " << meanSize[i] << std::endl;
    }
    out.close();
  }

  if (writeMembership) {
    uint64_t dataPrinted = 0;
    for (int p = 0; p < mpi_size; p++) {
      if (mpi_rank == p) {
        std::stringstream ss;
        ss << options.inputDir << "/membership_" << variant << options.outputSuffix;
        if (options.numRuns != 1)
          ss << "_" << options.currentRun;
        std::ofstream out(ss.str(), std::fstream::app);
        if (!out.is_open())
          throw std::runtime_error("Could not open output file for membership.");
      
        for (uint64_t i = 0; i < options.numLocalDataPoints; i++) {
          out << dataPrinted + i << " " << getMean(belongsToMean, i) << std::endl;
        }
        out.close();
      }
      else if (mpi_rank > p) {
        uint64_t fromp,top;
        std::tie(fromp,top) = mpi_range(options.numDataPoints, p);
        uint64_t beingPrinted = top - fromp;
        dataPrinted += beingPrinted;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  double end_time = MPI_Wtime();
  if (mpi_rank == 0)
    std::cout << "EXP " << options.currentRun << ": writeTime "
              << end_time - start_time << " seconds" << std::endl;
}


/*
 * Common parts of the k-Means algorithm.
 *
 * These implementations are generic in nature and are independent of
 * the two defined data structures (localized and non-localized).
 * When compiling a particular variant, support for the corresponding
 * data structure is handled by inlining the accessor functions defined
 * above into these templates.
 */

static inline void divideMeans(const struct Options &options,
                               double **meanValues,
                               uint64_t *meanSize)
{
  for (int i = 0; i < options.numMeans; i++) {
    for (int d = 0; d < options.dataDim; d++) {
      meanValues[i][d] = meanValues[i][d] / meanSize[i];
    }
  }
}

template<typename D, typename B>
static inline void initializeMeans(const struct Options &options,
                                   D *data,
                                   uint64_t *meanSize,
                                   uint64_t *meanSizeBuff,
                                   double **meanValues,
                                   double *meanValuesBuff,
                                   B *belongsToMean)
{
  std::fill(meanSizeBuff, meanSizeBuff + options.numMeans, 0);
  std::fill(meanValuesBuff, meanValuesBuff + options.numMeans * options.dataDim, 0.0);
  
  //ensure there is at least 1 data point at each mean
  int start = 0;
  if (mpi_rank == 0) {
    for (int i = 0; i < options.numMeans; i++) {
      setMean(belongsToMean, i, i);
      meanSizeBuff[i]++;
      
      for (int d = 0; d < options.dataDim; d++) {
        meanValuesBuff[i * options.dataDim + d] += getDataPoint(data, i, d);
      }
    }
    start = options.numMeans;
  }
  for (uint64_t i = start; i < options.numLocalDataPoints; i++) {
    int mean = rand() % options.numMeans;
    setMean(belongsToMean, i, mean);
    meanSizeBuff[mean]++;
    
    for (int d = 0; d < options.dataDim; d++) {
      meanValuesBuff[mean * options.dataDim + d] += getDataPoint(data, i, d);
    }
  }

  //communicate the inital values for the means and sizes of the means
  MPI_Allreduce(meanValuesBuff, meanValues[0], options.numMeans * options.dataDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(meanSizeBuff, meanSize, options.numMeans, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
}

inline std::vector<double> split_line(const std::string line) {
  std::vector<double> v;
  std::stringstream ss(line);
  std::string s;

  getline(ss, s, ' '); // skip point index
  while (getline(ss, s, ' ')) {
    v.push_back(std::stod(s));
  }

  return v;
}

inline void vector_to_mean_values(std::vector<std::vector<double>> v, double **meanValues) {
  for (unsigned int i = 0; i < v.size(); i++) {
    for (unsigned int j = 0; j < v[i].size(); j++) {
      meanValues[i][j] = v[i][j];
    }
  }
}

template<typename D, typename B>
inline void assign_initial_points(const struct Options &options, D *data, uint64_t *meanSize, uint64_t *meanSizeBuff, double **meanValues, double *meanValuesBuff, B *belongsToMean) {
  // TODO init: create other assign pattern (such as sequential assignment of points to clusters)
  int mean;

  for (uint64_t i = 0; i < options.numLocalDataPoints; i++) {
    mean = rand() % options.numMeans;
    setMean(belongsToMean, i, mean);
    meanSizeBuff[mean]++;
    
    for (int d = 0; d < options.dataDim; d++) {
      meanValuesBuff[mean * options.dataDim + d] += getDataPoint(data, i, d);
    }
  }

  //communicate the inital values for the means and sizes of the means
  MPI_Allreduce(meanValuesBuff, meanValues[0], options.numMeans * options.dataDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(meanSizeBuff, meanSize, options.numMeans, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
}

template<typename D, typename B>
static inline void initialize_means_from_file(const struct Options &options,
                                              D *data,
                                              uint64_t *meanSize,
                                              uint64_t *meanSizeBuff,
                                              double **meanValues,
                                              double *meanValuesBuff,
                                              B *belongsToMean) {
  // TODO INIT: check array contents
  std::fill(meanSizeBuff, meanSizeBuff + options.numMeans, 1);
  std::fill(meanValuesBuff, meanValuesBuff + options.numMeans * options.dataDim, 0.0);

  std::vector<std::vector<double>> v;
  std::string filename = options.inputDir;
  filename += + "/initial_means_" + std::to_string(options.meansSet) + ".txt";
  std::ifstream f(filename);
  std::string line;

  if (!f.is_open()) {
    std::cerr << "ERROR: unable to open file " << filename << std::endl;
  }
  if (mpi_rank == 0) {
    std::cout << "EXP " << options.currentRun << ": using initial means from file '"
              << filename << "'" << std::endl;
  }

  while (f >> line) {
    v.push_back(split_line(line));
  }

  vector_to_mean_values(v, meanValues);
  assign_initial_points(options, data, meanSize, meanSizeBuff, meanValues, meanValuesBuff, belongsToMean);
}

template<typename D, typename B>
static inline void recalcMeans(const struct Options &options,
                               D *data,
                               uint64_t *meanSize,
                               uint64_t *meanSizeBuff,
                               double **meanValues,
                               double *meanValuesBuff,
                               B *belongsToMean)
{
  std::fill(meanSizeBuff, meanSizeBuff + options.numMeans, 0);
  std::fill(meanValuesBuff,
            meanValuesBuff + options.numMeans * options.dataDim, 0.0);

  for (uint64_t i = 0; i < options.numLocalDataPoints; i++) {
    meanSizeBuff[getMean(belongsToMean, i)]++;

    for (int d = 0; d < options.dataDim; d++) {
      meanValuesBuff[getMean(belongsToMean, i) * options.dataDim + d] += getDataPoint(data, i, d);
    }
  }

  MPI_Allreduce(meanValuesBuff,meanValues[0], options.numMeans * options.dataDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(meanSizeBuff, meanSize, options.numMeans, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
}

static inline void updateMeans(const struct Options &options,
                               uint64_t *meanSize,
                               uint64_t *meanSizeBuff,
                               const uint64_t *oldMeanSize,
                               double **meanValues,
                               double *meanValuesBuff,
                               const double **oldMeanValuesSum)
{
  /* 
   * Let V_i0 and S_i0 be the value and size of mean i at time 0 (before this iteration)
   * Let V_ij0 and S_ij0 be the value and size of mean i at process j after reassigning the data
   * Then:
   *   S_i1 = sum S_ij0 - (p-1) * S_i0
   *   V_i1 = (sum (S_ij0 * V_ij0) - (p-1) * V_i0 * S_i0) / S_i1
   * where p is the total number of processes.
   *
   * Note: oldMeanValuesSum is filled with values of V_i0 * S_i0 (instead of just V_i0)
   */
  
  //communicate the sum of all current sizes and all current values multiplied by the size
  for (int mean = 0; mean < options.numMeans; mean++) {
    for (int d = 0; d < options.dataDim; d++) {
      meanValues[mean][d] = meanValues[mean][d] * meanSize[mean];
    }
  }
  MPI_Allreduce(meanValues[0], meanValuesBuff,
                options.numMeans * options.dataDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(meanSize, meanSizeBuff, options.numMeans, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  
  for (int mean = 0; mean < options.numMeans; mean++) {
    meanSize[mean] = meanSizeBuff[mean] - (mpi_size - 1) * oldMeanSize[mean];
    
    for (int d = 0; d < options.dataDim; d++) {
      meanValues[mean][d] = meanValuesBuff[mean * options.dataDim + d] - (mpi_size - 1) * oldMeanValuesSum[mean][d];
    }//Note: division by size will be done later
  }
}

template<typename D, typename B>
static inline void ensureNoEmptyClusters(const struct Options &options,
                                         uint64_t &specialReassigned,
                                         D *data,
                                         uint64_t *meanSize,
                                         double **meanValues,
                                         B *belongsToMean)
{
  //TODO: possibly parallelize this instead of letting process 0 do all the work? Is this even beneficial since the number of means should be fairly low so work needed is very limited.
  if (mpi_rank == 0) {
    uint64_t candidate = 0;
    specialReassigned = 0;
      
    for (int m = 0; m < options.numMeans; m++) {
      if (meanSize[m] == 0) {
        specialReassigned++;
          
        while (candidate < options.numLocalDataPoints) {
          if (meanSize[getMean(belongsToMean, candidate)] > 1) {
            meanSize[m] = 1;
            meanSize[getMean(belongsToMean, candidate)] -= 1;
            for (int d = 0; d < options.dataDim; d++) {
              meanValues[m][d] = getDataPoint(data, candidate, d);
              meanValues[getMean(belongsToMean, candidate)][d] -= getDataPoint(data, candidate, d);
            }
            setMean(belongsToMean, candidate, m);
            break;
          }
          candidate++;
        }
      }
    }
  }
}

static inline void recordOldMeans(const struct Options &options,
                                  const double **meanValues,
                                  double **oldMeanValues,
                                  const uint64_t *meanSize,
                                  uint64_t *oldMeanSize)
{
  //remember the old values and sizes of the means for checking the convergence Delta
  std::memcpy(oldMeanValues[0], meanValues[0],
         options.numMeans * options.dataDim * sizeof(double));
  std::memcpy(oldMeanSize, meanSize, options.numMeans * sizeof(uint64_t));
}

static inline bool meanHasMoved(const struct Options &options,
                                double **meanValues,
                                double **oldMeanValues)
{
  //now check to see if means have not moved:
  int iter = 0;
  for (iter = 0; iter < options.numMeans; iter++) {
    if (calculateDistance(meanValues[iter], oldMeanValues[iter], options.dataDim) >= options.convergenceDelta) 
      break;
  }
  if (iter == options.numMeans) {
    //no mean moved more than the convergence delta so we are done
    if (mpi_rank == 0) 
      std::cout << "Stopping iterations due to convergence delta." << std::endl;
    
    return false;
  }

  return true;
}

static inline void broadcastMeans(const struct Options &options,
                                  double **meanValues,
                                  uint64_t *meanSize)
{
  MPI_Bcast(meanValues[0], options.numMeans * options.dataDim,
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
  /* FIXME: shouldn't MPI_INT be the type corresponding to uint64_t? */
  MPI_Bcast(meanSize, options.numMeans, MPI_INT, 0, MPI_COMM_WORLD);
}

#endif /* __COMMON_H__ */
