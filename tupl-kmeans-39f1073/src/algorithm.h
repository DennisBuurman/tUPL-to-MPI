/*
 * MPI K-means implementation
 *
 * Author: Anne Hommelberg
 *
 * Code refactoring done by Kristian Rietveld, Leiden University.
 * 
 * Altered by Dennis Buurman to include new variants and functions
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

  // uint32_t count = 0;
  // if (mpi_rank == 0) {
  //   std::cout << "LOCAL POINTS: " <<  numLocalDataPoints;
  // }

  double oldDistance, newDistance;
  for (uint64_t x = 0; x < numLocalDataPoints; x++) {
    oldDistance = calculateDistance(meanValues[getMean(belongsToMean, x)], getData(data, x), dataDim);
    for (int m = 0; m < numMeans; m++) {
      if (m != getMean(belongsToMean, x)) {
        newDistance = calculateDistance(meanValues[m], getData(data, x), dataDim);
        if (newDistance < oldDistance) {
          localReassigned++;
          // count++;
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

  // if (mpi_rank == 0) {
  //   std::cout << "; REASSIGNED: " << count << std::endl;
  // }
}

template <typename D, typename B>
void kmeansRecalc(struct Options &options, const std::string &variant,
                  D *data, B *belongsToMean)
{
  double start_time, end_time;

  // Time the individual parts of computation
  double init_time, it_start_time, it_reassign_time, it_update_time, it_end_time;
  double cum_reassign_time = 0.0, cum_update_time = 0.0, cum_it_end_time = 0.0;

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

  // Initialize set of means from file if provided
  if (options.meansFlag) {
    initialize_means_from_file(options, data, meanSize, meanSizeBuff,
                               meanValues, meanValuesBuff, belongsToMean);
  } else {
    initializeMeans(options, data, meanSize, meanSizeBuff,
                    meanValues, meanValuesBuff, belongsToMean);
    divideMeans(options, meanValues, meanSize);
  }

  recordOldMeans(options,
                 const_cast<const double **>(meanValues), oldMeanValues,
                 const_cast<const uint64_t *>(meanSize), oldMeanSize);

  init_time = MPI_Wtime(); // time from start to end of initialization

  //execute the algorithm (this is where the magic happens)
  //NOTE: this is only one of several possible ways to code this
  uint64_t reassigned = threshold + 1;
  uint64_t localReassigned, specialReassigned = 0;
  int cycles = 0;
  while (reassigned > threshold) {
    it_start_time = MPI_Wtime();  // time start of each iteration

    cycles++;
    localReassigned = 0;

    reassignLocalDataPoints(options, localReassigned, data,
                            meanSize, meanValues, belongsToMean);
    
    it_reassign_time = MPI_Wtime();  // time after point reassignment

    // now recalculate means (by communication) and check whether we are done
    recalcMeans(options, data, meanSize, meanSizeBuff,
                meanValues, meanValuesBuff, belongsToMean); 
    
    it_update_time = MPI_Wtime();  // time after updating the means (buffering+communication)

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

    it_end_time = MPI_Wtime();  // time end of each iteration (mean division)

    // Update timing variables
    cum_reassign_time += it_reassign_time - it_start_time;
    cum_update_time += it_update_time - it_reassign_time;
    cum_it_end_time += it_end_time - it_update_time;
  }
  
  //log the end time and elapsed time
  end_time = MPI_Wtime();
  std::cout << "EXP " << options.currentRun << ": " << printPreamble()
            << " calculationTime "
            << end_time - start_time << " seconds" << std::endl;
  
  //log the time spent in each component
  // if (mpi_rank == 0) {
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " initTime " << init_time - start_time << " s" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " calculationTime " << end_time - start_time << " s | "
              << (end_time - start_time) / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " reassignTime " << cum_reassign_time << " s | "
              << cum_reassign_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " communicationTime "  << cum_update_time << " s | "
              << cum_update_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " meanRecalculationTime " << cum_it_end_time << " s | "
              << cum_it_end_time / cycles << " s averaged" << std::endl;
  // }

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

  // Time the individual parts of computation
  double init_time, it_start_time, it_reassign_time, it_update_time, it_end_time;
  double cum_reassign_time = 0.0, cum_update_time = 0.0, cum_it_end_time = 0.0;

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

  // Initialize set of means from file if provided
  if (options.meansFlag) {
    initialize_means_from_file(options, data, meanSize, meanSizeBuff,
                               meanValues, meanValuesBuff, belongsToMean);
  } else {
    initializeMeans(options, data, meanSize, meanSizeBuff,
                    meanValues, meanValuesBuff, belongsToMean);
    divideMeans(options, meanValues, meanSize);
  }

  //remember the old values and sizes of the means for recalculation
  memcpy(oldMeanValuesSum[0], meanValues[0], numMeans * dataDim * sizeof(double));

  recordOldMeans(options,
                 const_cast<const double **>(meanValues), oldMeanValues,
                 const_cast<const uint64_t *>(meanSize), oldMeanSize);
  
  init_time = MPI_Wtime(); // time from start to end of initialization

  //execute the algorithm (this is where the magic happens)
  //NOTE: this is only one of several possible ways to code this
  
  uint64_t reassigned = threshold + 1;
  uint64_t localReassigned, specialReassigned = 0;
  uint64_t reassignedSinceRecalculation = 0;
  int cycles = 0;
  while (reassigned > threshold) {
    it_start_time = MPI_Wtime();  // time start of each iteration

    cycles++;
    localReassigned = 0;

    reassignLocalDataPoints(options, localReassigned, data,
                            meanSize, meanValues, belongsToMean);
    
    it_reassign_time = MPI_Wtime();  // time after reassignment

    // communicate whether all processes have converged
    MPI_Allreduce(&localReassigned, &reassigned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //TODO: usefull to add a break if no reassignments are made?
    
    // now communicate the new values for the means between processes
    reassignedSinceRecalculation += reassigned;
    if (reassignedSinceRecalculation > REASSIGN_LIMIT || cycles < 2) {
      // initializing from file requires one recalculation before means can be updated.
      // this is because of oldMeanSize and oldMeanValuesSum
      // TODO: only do cycles check when mean set flag is set!
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
    
    it_update_time = MPI_Wtime();  // time after updating the means (buffering+communication)

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
    
    it_end_time = MPI_Wtime();  // time at end of iteration (mean division)

    // Update timing variables
    cum_reassign_time += it_reassign_time - it_start_time;
    cum_update_time += it_update_time - it_reassign_time;
    cum_it_end_time += it_end_time - it_update_time;
  }
  
  //log the end time and elapsed time
  end_time = MPI_Wtime();
  std::cout << "EXP " << options.currentRun << ": " << printPreamble()
            << " calculationTime "
            << end_time - start_time << " seconds" << std::endl;
  
  //log the time spent in each component
  // if (mpi_rank == 0) {
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " initTime " << init_time - start_time << " s" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " calculationTime " << end_time - start_time << " s | "
              << (end_time - start_time) / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " reassignTime " << cum_reassign_time << " s | "
              << cum_reassign_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " communicationTime "  << cum_update_time << " s | "
              << cum_update_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " meanRecalculationTime " << cum_it_end_time << " s | "
              << cum_it_end_time / cycles << " s averaged" << std::endl;
  // }

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

/* Additions to file */

// Initialize the local partitions of meanValues and meanSize.
template <typename D, typename B>
static inline void initLocal(struct Options &options,
                             D *data,
                             B *belongsToMean,
                             uint64_t * meanSizeLocal,
                             double **meanValuesLocal)
{
  const int numMeans(options.numMeans);
  const uint64_t numLocalDataPoints(options.numLocalDataPoints);
  const int dataDim(options.dataDim);

  std::fill(meanSizeLocal, meanSizeLocal + numMeans, 0);
  std::fill(&meanValuesLocal[0][0], &meanValuesLocal[0][0] + sizeof(meanValuesLocal), 0.0);

  // Add all local points to their corresponding mean and size arrays
  int mean;
  for (uint64_t x = 0; x < numLocalDataPoints; x++) {
    mean = getMean(belongsToMean, x);
    meanSizeLocal[mean]++;
    for (int d = 0; d < dataDim; d++) {
      meanValuesLocal[mean][d] += getDataPoint(data, x, d);
    }
  }

  // Perform mean division to receive initial local mean values
  for (int mean = 0; mean < numMeans; mean++) {
    for (int d = 0; d < dataDim; d++) {
      meanValuesLocal[mean][d] = meanValuesLocal[mean][d] / meanSizeLocal[mean];
    }
  }
}

// Initialize the local partition of mean values only.
// Local partition of mean values contains sum of all local point coords instead of mean.
template <typename D, typename B>
static inline void initLocalValues(struct Options &options,
                             D *data,
                             B *belongsToMean,
                             double **meanValuesLocal)
{
  const int numMeans(options.numMeans);
  const uint64_t numLocalDataPoints(options.numLocalDataPoints);
  const int dataDim(options.dataDim);

  std::fill(&meanValuesLocal[0][0], &meanValuesLocal[0][0] + sizeof(meanValuesLocal), 0.0);

  // Sum coordinates of all local points in the correct mean
  int mean;
  for (uint64_t x = 0; x < numLocalDataPoints; x++) {
    mean = getMean(belongsToMean, x);
    for (int d = 0; d < dataDim; d++) {
      meanValuesLocal[mean][d] += getDataPoint(data, x, d);
    }
  }
}

// Reassign but also keep track of mean values and mean sizes for local data points
template <typename D, typename B>
static inline void reassignLocalDataPointsMlevel(struct Options &options,
                                           uint64_t &localReassigned,
                                           D *data,
                                           uint64_t *meanSize,
                                           uint64_t * meanSizeLocal,
                                           double **meanValues,
                                           double **meanValuesLocal,
                                           B *belongsToMean)
{
  const int numMeans(options.numMeans);
  const uint64_t numLocalDataPoints(options.numLocalDataPoints);
  const int dataDim(options.dataDim);

  double oldDistance, newDistance;

  // uint32_t count = 0;
  // if (mpi_rank == 0) {
  //   std::cout << "LOCAL POINTS: " <<  numLocalDataPoints;
  // }

  for (uint64_t x = 0; x < numLocalDataPoints; x++) {
    oldDistance = calculateDistance(meanValues[getMean(belongsToMean, x)], getData(data, x), dataDim);
    for (int m = 0; m < numMeans; m++) {
      if (m != getMean(belongsToMean, x)) {
        newDistance = calculateDistance(meanValues[m], getData(data, x), dataDim);
        if (newDistance < oldDistance) {
          localReassigned++;
          // count++;
          const int oldMean(getMean(belongsToMean, x));
          for (int d = 0; d < dataDim; d++) {
            /* Updating the local copy of global mean values and sizes affects point reassignments, but has little effect on end results! */
            // meanValues[oldMean][d] = (meanValues[oldMean][d] * meanSize[oldMean] - getDataPoint(data, x, d)) / (meanSize[oldMean] - 1);
            // meanValues[m][d] = (meanValues[m][d] * meanSize[m] + getDataPoint(data, x, d)) / (meanSize[m] + 1);

            meanValuesLocal[oldMean][d] = (meanValuesLocal[oldMean][d] * meanSizeLocal[oldMean] - getDataPoint(data, x, d)) / (meanSizeLocal[oldMean] - 1);
            meanValuesLocal[m][d] = (meanValuesLocal[m][d] * meanSizeLocal[m] + getDataPoint(data, x, d)) / (meanSizeLocal[m] + 1);
          }
          /* Updating the local copy of global mean values and sizes affects point reassignments, but has little effect on end results! */
          // meanSize[oldMean] -= 1;
          // meanSize[m] += 1;
          setMean(belongsToMean, x, m);
          meanSizeLocal[oldMean] -= 1;
          meanSizeLocal[m] += 1;
          oldDistance = newDistance; 
        }
      }
    }
  }

  // if (mpi_rank == 0) {
  //   std::cout << "; REASSIGNED: " << count << std::endl;
  // }
}

// Original recalc variant without updating local version of global means
// Removing the updates locally affects point reassignment, but has little effect on end results
template <typename D, typename B>
static inline void reassignLocalDataPointsNoUpdates(struct Options &options,
                                           uint64_t &localReassigned,
                                           D *data,
                                           uint64_t *meanSize,
                                           double **meanValues,
                                           B *belongsToMean)
{
  const int numMeans(options.numMeans);
  const uint64_t numLocalDataPoints(options.numLocalDataPoints);
  const int dataDim(options.dataDim);

  // uint32_t count = 0;
  // if (mpi_rank == 0) {
  //   std::cout << "LOCAL POINTS: " <<  numLocalDataPoints;
  // }

  double oldDistance, newDistance;
  for (uint64_t x = 0; x < numLocalDataPoints; x++) {
    oldDistance = calculateDistance(meanValues[getMean(belongsToMean, x)], getData(data, x), dataDim);
    for (int m = 0; m < numMeans; m++) {
      if (m != getMean(belongsToMean, x)) {
        newDistance = calculateDistance(meanValues[m], getData(data, x), dataDim);
        if (newDistance < oldDistance) {
          localReassigned++;
          // count++;
          setMean(belongsToMean, x, m);
          oldDistance = newDistance; 
        }
      }
    }
  }

  // if (mpi_rank == 0) {
  //   std::cout << "; REASSIGNED: " << count << std::endl;
  // }
}

// Ignore write-dependency of meanValues (meanSize dependency)
// reassign each of the local data points, updating the mean values as we go
// meanSize is kept strongly consistent. Each update is directly synchronized.
// Alternative version of reassignLocalDataPoints function.
template <typename D, typename B>
static inline void reassignSyncSize(struct Options &options,
                                    uint64_t &localReassigned,
                                    D *data,
                                    uint64_t *meanSize,
                                    double **meanValues,
                                    double **meanValuesLocal, // local partition
                                    B *belongsToMean)
{
  const int numMeans(options.numMeans);
  const uint64_t numLocalDataPoints(options.numLocalDataPoints);
  const int dataDim(options.dataDim);

  int64_t *meanSizeChange = new int64_t[numMeans];

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
            meanValuesLocal[oldMean][d] -= getDataPoint(data, x, d);
            meanValuesLocal[m][d] += getDataPoint(data, x, d);
          }
          meanSize[oldMean] -= 1;
          meanSize[m] += 1;

          /* Communicate change and update meanSize */
          std::fill(meanSizeChange, meanSizeChange + numMeans, 0);

          meanSizeChange[oldMean] -= 1;
          meanSizeChange[m] += 1;

          MPI_Allreduce(MPI_IN_PLACE, meanSizeChange, options.numMeans, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

          for (int i = 0; i < options.numMeans; i++) {
            meanSize[i] += meanSizeChange[i];
          }

          setMean(belongsToMean, x, m);
          oldDistance = newDistance; 
        }
      }
    }
  }

  delete[] meanSizeChange;
  meanSizeChange = NULL;
}

// Recalculation variant that ignores meanSize dependency from meanValues
// meanSize is kept strongly consistent, while meanValues is updated after each iteration
template <typename D, typename B>
void kmeansRecalcNoDependencies(struct Options &options, const std::string &variant,
                  D *data, B *belongsToMean)
{
  double start_time, end_time;

  // Time the individual parts of computation
  double init_time, it_start_time, it_reassign_time, it_update_time, it_end_time;
  double cum_reassign_time = 0.0, cum_update_time = 0.0, cum_it_end_time = 0.0;

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
  double ** meanValuesLocal = allocate2dDoubleArray(numMeans,dataDim); // local partition
  double * meanValuesBuff = new double[numMeans * dataDim];
  double threshold = numDataPoints * options.thresholdMultiplier;

  if (options.meansFlag) {
    initialize_means_from_file(options, data, meanSize, meanSizeBuff,
                               meanValues, meanValuesBuff, belongsToMean);
  } else {
    initializeMeans(options, data, meanSize, meanSizeBuff,
                    meanValues, meanValuesBuff, belongsToMean);
    divideMeans(options, meanValues, meanSize);
  }

  recordOldMeans(options,
                 const_cast<const double **>(meanValues), oldMeanValues,
                 const_cast<const uint64_t *>(meanSize), oldMeanSize);

  initLocalValues(options, data, belongsToMean, meanValuesLocal); // NEW: init local partition

  init_time = MPI_Wtime();

  //execute the algorithm (this is where the magic happens)
  //NOTE: this is only one of several possible ways to code this
  uint64_t reassigned = threshold + 1;
  uint64_t localReassigned, specialReassigned = 0;
  int cycles = 0;
  while (reassigned > threshold) {
    it_start_time = MPI_Wtime();  // at start of iteration

    cycles++;
    localReassigned = 0;

    // NEW: COMMUNICATE SIZE AFTER EACH REASSIGNMENT
    reassignSyncSize(options, localReassigned, data, meanSize,
                     meanValues, meanValuesLocal, belongsToMean);
    
    it_reassign_time = MPI_Wtime();  // after reassignment

    // NEW: ONLY COMMUNICATE meanValuesBuff
    // now recalculate means (by communication) and check whether we are done
    recalcMeansValuesOnly(options, data, meanValues, meanValuesBuff, meanValuesLocal, belongsToMean);
    
    it_update_time = MPI_Wtime();  // after updating the means

    // NEW: disabled (TODO: check if needed)
    // //check to avoid empty clusters, if so, reassign a point to fill them
    // ensureNoEmptyClusters(options, specialReassigned,
    //                       data, meanSize, meanValues, belongsToMean); 
    // broadcastMeans(options, meanValues, meanSize);

    divideMeans(options, meanValues, meanSize);
    
    if (!meanHasMoved(options, meanValues, oldMeanValues))
      break;
    
    recordOldMeans(options,
                   const_cast<const double **>(meanValues), oldMeanValues,
                   const_cast<const uint64_t *>(meanSize), oldMeanSize);

    // communicate whether all processes have converged
    MPI_Allreduce(&localReassigned, &reassigned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    it_end_time = MPI_Wtime();  // at end of iteration

    // Update timing variables
    cum_reassign_time += it_reassign_time - it_start_time;
    cum_update_time += it_update_time - it_reassign_time;
    cum_it_end_time += it_end_time - it_update_time;
  }
  
  //log the end time and elapsed time
  end_time = MPI_Wtime();
  std::cout << "EXP " << options.currentRun << ": " << printPreamble()
            << " calculationTime "
            << end_time - start_time << " seconds" << std::endl;
  
  //log the time spent in each component
  // if (mpi_rank == 0) {
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " initTime " << init_time - start_time << " s" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " calculationTime " << end_time - start_time << " s | "
              << (end_time - start_time) / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " reassignTime " << cum_reassign_time << " s | "
              << cum_reassign_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " communicationTime "  << cum_update_time << " s | "
              << cum_update_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " meanRecalculationTime " << cum_it_end_time << " s | "
              << cum_it_end_time / cycles << " s averaged" << std::endl;
  // }

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
  delete[] meanValuesLocal[0]; // local partition
  delete[] meanValuesLocal; // local partition
}

// Recalculation variant that adds and tracks local partitions.
// This way, we can communicate size*mean instead of sum(localdatapoints)
template <typename D, typename B>
void kmeansRecalcMlevel(struct Options &options, const std::string &variant,
                  D *data, B *belongsToMean)
{
  double start_time, end_time;

  // Time the individual parts of computation
  double init_time, it_start_time, it_reassign_time, it_update_time, it_end_time;
  double cum_reassign_time = 0.0, cum_update_time = 0.0, cum_it_end_time = 0.0;

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
  uint64_t * meanSizeLocal = new uint64_t[numMeans]; // NEW
  double ** meanValues = allocate2dDoubleArray(numMeans,dataDim);
  double ** oldMeanValues = allocate2dDoubleArray(numMeans,dataDim);
  double ** meanValuesLocal = allocate2dDoubleArray(numMeans,dataDim); // NEW
  double * meanValuesBuff = new double[numMeans * dataDim];
  double threshold = numDataPoints * options.thresholdMultiplier;

  if (options.meansFlag) {
    initialize_means_from_file(options, data, meanSize, meanSizeBuff,
                               meanValues, meanValuesBuff, belongsToMean);
  } else {
    initializeMeans(options, data, meanSize, meanSizeBuff,
                    meanValues, meanValuesBuff, belongsToMean);
    divideMeans(options, meanValues, meanSize);
  }


  recordOldMeans(options,
                 const_cast<const double **>(meanValues), oldMeanValues,
                 const_cast<const uint64_t *>(meanSize), oldMeanSize);

  initLocal(options, data, belongsToMean, meanSizeLocal, meanValuesLocal); // NEW

  init_time = MPI_Wtime();

  //execute the algorithm (this is where the magic happens)
  //NOTE: this is only one of several possible ways to code this
  uint64_t reassigned = threshold + 1;
  uint64_t localReassigned, specialReassigned = 0;
  const uint64_t REASSIGN_LIMIT = 1000000000; // NEW
  uint64_t reassignedSinceRecalculation = 0; // NEW
  int cycles = 0;
  while (reassigned > threshold) {
    it_start_time = MPI_Wtime();  // at start of iteration

    cycles++;
    localReassigned = 0;

    reassignLocalDataPointsMlevel(options, localReassigned, data, meanSize, meanSizeLocal,
                                  meanValues, meanValuesLocal, belongsToMean); // NEW
    
    it_reassign_time = MPI_Wtime();  // after reassignment

    // now recalculate means (by communication) and check whether we are done
    MPI_Allreduce(&localReassigned, &reassigned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // NEW
    reassignedSinceRecalculation += reassigned;
    if (reassignedSinceRecalculation > REASSIGN_LIMIT) {
      if (mpi_rank == 0) {
        std::cout << "using i-level due to rounding errors" << std::endl;
      }
      //we've reassigned enough points to risk significant rounding errors, so do a full recalculation
      reassignedSinceRecalculation = 0;
      
      recalcMeans(options, data, meanSize, meanSizeBuff,
                  meanValues, meanValuesBuff, belongsToMean);
    }
    else {
      recalcMeansMlevel(options, data, meanSize, meanSizeBuff, meanSizeLocal,
                      meanValues, meanValuesBuff, meanValuesLocal, belongsToMean); // NEW
    }
    
    it_update_time = MPI_Wtime();  // after updating the means

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

    it_end_time = MPI_Wtime();  // at end of iteration

    // Update timing variables
    cum_reassign_time += it_reassign_time - it_start_time;
    cum_update_time += it_update_time - it_reassign_time;
    cum_it_end_time += it_end_time - it_update_time;
  }
  
  //log the end time and elapsed time
  end_time = MPI_Wtime();
  std::cout << "EXP " << options.currentRun << ": " << printPreamble()
            << " calculationTime "
            << end_time - start_time << " seconds" << std::endl;
  
  //log the time spent in each component
  // if (mpi_rank == 0) {
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " initTime " << init_time - start_time << " s" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " calculationTime " << end_time - start_time << " s | "
              << (end_time - start_time) / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " reassignTime " << cum_reassign_time << " s | "
              << cum_reassign_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " communicationTime "  << cum_update_time << " s | "
              << cum_update_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " meanRecalculationTime " << cum_it_end_time << " s | "
              << cum_it_end_time / cycles << " s averaged" << std::endl;
  // }

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
  delete[] meanSizeLocal; // NEW
  delete[] meanValuesLocal[0]; // NEW
  delete[] meanValuesLocal; // NEW
}

// Recalculation variant that removes updating the local copy of global meanValues and meanSize arrays during point reassignment.
template <typename D, typename B>
void kmeansRecalcNoUpdates(struct Options &options, const std::string &variant,
                  D *data, B *belongsToMean)
{
  double start_time, end_time;

  // Time the individual parts of computation
  double init_time, it_start_time, it_reassign_time, it_update_time, it_end_time;
  double cum_reassign_time = 0.0, cum_update_time = 0.0, cum_it_end_time = 0.0;

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

  if (options.meansFlag) {
    initialize_means_from_file(options, data, meanSize, meanSizeBuff,
                               meanValues, meanValuesBuff, belongsToMean);
  } else {
    initializeMeans(options, data, meanSize, meanSizeBuff,
                    meanValues, meanValuesBuff, belongsToMean);
    divideMeans(options, meanValues, meanSize);
  }


  recordOldMeans(options,
                 const_cast<const double **>(meanValues), oldMeanValues,
                 const_cast<const uint64_t *>(meanSize), oldMeanSize);

  init_time = MPI_Wtime();

  //execute the algorithm (this is where the magic happens)
  //NOTE: this is only one of several possible ways to code this
  uint64_t reassigned = threshold + 1;
  uint64_t localReassigned, specialReassigned = 0;
  int cycles = 0;
  while (reassigned > threshold) {
    it_start_time = MPI_Wtime();  // at start of iteration

    cycles++;
    localReassigned = 0;

    // IM version (no updates on local copy of global means)
    reassignLocalDataPointsNoUpdates(options, localReassigned, data,
                            meanSize, meanValues, belongsToMean);
    
    it_reassign_time = MPI_Wtime();  // after reassignment

    // now recalculate means (by communication) and check whether we are done
    recalcMeans(options, data, meanSize, meanSizeBuff,
                meanValues, meanValuesBuff, belongsToMean); 
    
    it_update_time = MPI_Wtime();  // after updating the means

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

    it_end_time = MPI_Wtime();  // at end of iteration

    // Update timing variables
    cum_reassign_time += it_reassign_time - it_start_time;
    cum_update_time += it_update_time - it_reassign_time;
    cum_it_end_time += it_end_time - it_update_time;
  }
  
  //log the end time and elapsed time
  end_time = MPI_Wtime();
  std::cout << "EXP " << options.currentRun << ": " << printPreamble()
            << " calculationTime "
            << end_time - start_time << " seconds" << std::endl;
  
  //log the time spent in each component
  // if (mpi_rank == 0) {
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " initTime " << init_time - start_time << " s" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " calculationTime " << end_time - start_time << " s | "
              << (end_time - start_time) / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " reassignTime " << cum_reassign_time << " s | "
              << cum_reassign_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " communicationTime "  << cum_update_time << " s | "
              << cum_update_time / cycles << " s averaged" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " meanRecalculationTime " << cum_it_end_time << " s | "
              << cum_it_end_time / cycles << " s averaged" << std::endl;
  // }

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