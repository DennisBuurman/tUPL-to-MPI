/*
 * MPI K-means implementation
 *
 * Author: Anne Hommelberg
 *
 * Code refactoring done by Kristian Rietveld, Leiden University.
 */

#include "mpi-utils.h"
#include "common.h"

#include <iostream>
#include <string>

// reassign each of the local data points, updating the mean values as we go
template <typename D, typename B>
static inline void reassignLocalDataPoints(struct Options &options,
                                           uint64_t &localReassigned,
                                           D *data,
                                           uint64_t *meanSize,
                                           double **meanValues,
                                           B *belongsToMean)
{
  const int numMeans(options.numMeans);
  const uint64_t numLocalDataPoints(options.numLocalDataPoints);
  const int dataDim(options.dataDim);

  double oldDistance, newDistance;
  for (uint64_t x = 0; x < numLocalDataPoints; x++) {
    oldDistance = calculateDistance(meanValues[getMean(belongsToMean, x)], getData(data, x), dataDim);
    for (int m = 0; m < numMeans; m++) {
      if (m != getMean(belongsToMean, x)) {
        newDistance = calculateDistance(meanValues[m], getData(data, x), dataDim);
        if (newDistance < oldDistance) {
          localReassigned++;
          const int oldMean(getMean(belongsToMean, x));
          for (int d = 0; d < dataDim; d++) {
            meanValues[oldMean][d] = (meanValues[oldMean][d] * meanSize[oldMean] - getDataPoint(data, x, d)) / (meanSize[oldMean] - 1);
            meanValues[m][d] = (meanValues[m][d] * meanSize[m] + getDataPoint(data, x, d)) / (meanSize[m] + 1);
          }
          meanSize[oldMean] -= 1;
          meanSize[m] += 1;
          setMean(belongsToMean, x, m);
          oldDistance = newDistance; 
        }
      }
    }
  }
}

template <typename D, typename B>
void kmeansRecalc(struct Options &options, const std::string &variant,
                  D *data, B *belongsToMean)
{
  double start_time, end_time;

  //define local variables for the often used options
  const int numMeans(options.numMeans);
  const uint64_t numDataPoints(options.numDataPoints);
  const int dataDim(options.dataDim);
  
  initRandom(options);

  //log the start time
  MPI_Barrier(MPI_COMM_WORLD);
  start_time = MPI_Wtime();
  
  //Initialization for algorithm
  uint64_t * meanSize = new uint64_t[numMeans];
  uint64_t * oldMeanSize = new uint64_t[numMeans];
  uint64_t * meanSizeBuff = new uint64_t[numMeans];
  double ** meanValues = allocate2dDoubleArray(numMeans,dataDim);
  double ** oldMeanValues = allocate2dDoubleArray(numMeans,dataDim);
  double * meanValuesBuff = new double[numMeans * dataDim];
  double threshold = numDataPoints * options.thresholdMultiplier;
  
  initializeMeans(options, data, meanSize, meanSizeBuff,
                  meanValues, meanValuesBuff, belongsToMean);
  divideMeans(options, meanValues, meanSize);

  recordOldMeans(options,
                 const_cast<const double **>(meanValues), oldMeanValues,
                 const_cast<const uint64_t *>(meanSize), oldMeanSize);
  
  //execute the algorithm (this is where the magic happens)
  //NOTE: this is only one of several possible ways to code this
  
  uint64_t reassigned = threshold + 1;
  uint64_t localReassigned, specialReassigned = 0;
  int cycles = 0;
  while (reassigned > threshold) {
    cycles++;
    localReassigned = 0;
    
    reassignLocalDataPoints(options, localReassigned, data,
                            meanSize, meanValues, belongsToMean);
    
    // now recalculate means (by communication) and check whether we are done
    recalcMeans(options, data, meanSize, meanSizeBuff,
                meanValues, meanValuesBuff, belongsToMean);
    
    //check to avoid empty clusters, if so, reassign a point to fill them
    ensureNoEmptyClusters(options, specialReassigned,
                          data, meanSize, meanValues, belongsToMean);

    broadcastMeans(options, meanValues, meanSize);
    divideMeans(options, meanValues, meanSize);
    
    if (!meanHasMoved(options, meanValues, oldMeanValues))
      break;
    
    recordOldMeans(options,
                   const_cast<const double **>(meanValues), oldMeanValues,
                   const_cast<const uint64_t *>(meanSize), oldMeanSize);

    // communicate whether all processes have converged
    MPI_Allreduce(&localReassigned, &reassigned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  
  //log the end time and elapsed time
  end_time = MPI_Wtime();
  std::cout << "EXP " << options.currentRun << ": " << printPreamble()
            << " calculationTime "
            << end_time - start_time << " seconds" << std::endl;
  
  if (mpi_rank == 0)
    std::cout << "EXP " << options.currentRun << ": iterations "
              << cycles << std::endl;

  writeData(options, variant, meanValues, meanSize, belongsToMean);

  delete[] meanSize;
  delete[] oldMeanSize;
  delete[] meanSizeBuff;
  delete[] meanValues[0];
  delete[] meanValues;
  delete[] oldMeanValues[0];
  delete[] oldMeanValues;
  delete[] meanValuesBuff;
}


const uint64_t REASSIGN_LIMIT = 1000000000;  //used for testing, should be set so no rounding errors occur, doubles should have approx 15 digit precision


template <typename D, typename B>
void kmeansIncremental(struct Options &options, const std::string &variant,
                       D *data, B *belongsToMean)
{
  double start_time, end_time;

  const int numMeans(options.numMeans);
  const uint64_t numDataPoints(options.numDataPoints);
  const int dataDim(options.dataDim);

  initRandom(options);
  
  //log the start time
  MPI_Barrier(MPI_COMM_WORLD);
  start_time = MPI_Wtime();
  
  //Initialization for algorithm
  uint64_t * meanSize = new uint64_t[numMeans];
  uint64_t * oldMeanSize = new uint64_t[numMeans];
  uint64_t * meanSizeBuff = new uint64_t[numMeans];
  double ** meanValues = allocate2dDoubleArray(numMeans,dataDim);
  double ** oldMeanValuesSum = allocate2dDoubleArray(numMeans,dataDim);
  double ** oldMeanValues = allocate2dDoubleArray(numMeans,dataDim);
  double * meanValuesBuff = new double[numMeans * dataDim];
  double threshold = numDataPoints * options.thresholdMultiplier;
  //TODO make meanValuesBuff a 2d array as well since this is just inconsistent
  
  initializeMeans(options, data, meanSize, meanSizeBuff,
                  meanValues, meanValuesBuff, belongsToMean);
  
  //remember the old values and sizes of the means for recalculation
  memcpy(oldMeanValuesSum[0], meanValues[0], numMeans * dataDim * sizeof(double));

  divideMeans(options, meanValues, meanSize);
  
  recordOldMeans(options,
                 const_cast<const double **>(meanValues), oldMeanValues,
                 const_cast<const uint64_t *>(meanSize), oldMeanSize);
  
  
  //execute the algorithm (this is where the magic happens)
  //NOTE: this is only one of several possible ways to code this
  
  uint64_t reassigned = threshold + 1;
  uint64_t localReassigned, specialReassigned = 0;
  uint64_t reassignedSinceRecalculation = 0;
  int cycles = 0;
  while (reassigned > threshold) {
    cycles++;
    localReassigned = 0;
    
    reassignLocalDataPoints(options, localReassigned, data,
                            meanSize, meanValues, belongsToMean);
    
    // communicate whether all processes have converged
    MPI_Allreduce(&localReassigned, &reassigned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //TODO: usefull to add a break if no reassignments are made?
    
    // now communicate the new values for the means between processes
    reassignedSinceRecalculation += reassigned;
    if (reassignedSinceRecalculation > REASSIGN_LIMIT) {
      if (mpi_rank == 0) {
        std::cout << "recalculating due to rounding errors" << std::endl;
      }
      //we've reassigned enough points to risk significant rounding errors, so do a full recalculation
      reassignedSinceRecalculation = 0;
      
      recalcMeans(options, data, meanSize, meanSizeBuff,
                  meanValues, meanValuesBuff, belongsToMean);
    }
    else {
      //no risk for rounding errors yet, so just recalculate using the previous values
      updateMeans(options,
                  meanSize, meanSizeBuff, const_cast<const uint64_t *>(oldMeanSize),
                  meanValues, meanValuesBuff, const_cast<const double **>(oldMeanValuesSum));
    }
    
    //check to avoid empty clusters, if so, reassign a point to fill them
    ensureNoEmptyClusters(options, specialReassigned,
                          data, meanSize, meanValues, belongsToMean);

    broadcastMeans(options, meanValues, meanSize);
    
    memcpy(oldMeanValuesSum[0], meanValues[0], numMeans * dataDim * sizeof(double));
    
    divideMeans(options, meanValues, meanSize);
    
    if (!meanHasMoved(options, meanValues, oldMeanValues))
      break;
    
    recordOldMeans(options,
                   const_cast<const double **>(meanValues), oldMeanValues,
                   const_cast<const uint64_t *>(meanSize), oldMeanSize);
  }
  
  //log the end time and elapsed time
  end_time = MPI_Wtime();
  std::cout << "EXP " << options.currentRun << ": " << printPreamble()
            << " calculationTime "
            << end_time - start_time << " seconds" << std::endl;
  
  if (mpi_rank == 0)
    std::cout << "EXP " << options.currentRun << ": iterations "
              << cycles << std::endl;

  writeData(options, variant, meanValues, meanSize, belongsToMean);

  delete[] meanSize;
  delete[] oldMeanSize;
  delete[] meanSizeBuff;
  delete[] meanValues[0];
  delete[] meanValues;
  delete[] oldMeanValuesSum[0];
  delete[] oldMeanValuesSum;
  delete[] oldMeanValues[0];
  delete[] oldMeanValues;
  delete[] meanValuesBuff;
}
