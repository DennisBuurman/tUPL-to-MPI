/*
 * MPI K-means implementation
 *
 * Author: Anne Hommelberg
 *
 * Code refactoring done by Kristian Rietveld, Leiden University.
 * 
 * Altered by Dennis Buurman to include new variants
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
          // // *** TURN OFF FOR IM variant
          // const int oldMean(getMean(belongsToMean, x));
          // for (int d = 0; d < dataDim; d++) {
          //   meanValues[oldMean][d] = (meanValues[oldMean][d] * meanSize[oldMean] - getDataPoint(data, x, d)) / (meanSize[oldMean] - 1);
          //   meanValues[m][d] = (meanValues[m][d] * meanSize[m] + getDataPoint(data, x, d)) / (meanSize[m] + 1);
          // }
          // meanSize[oldMean] -= 1;
          // meanSize[m] += 1;
          // // ***
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
  uint64_t * meanSizeLocal = new uint64_t[numMeans]; // M level recalc
  double ** meanValues = allocate2dDoubleArray(numMeans,dataDim);
  double ** oldMeanValues = allocate2dDoubleArray(numMeans,dataDim);
  double ** meanValuesLocal = allocate2dDoubleArray(numMeans,dataDim); // M level recalc
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

  initLocal(options, data, belongsToMean, meanSizeLocal, meanValuesLocal); // M level recalc

  init_time = MPI_Wtime();

  //execute the algorithm (this is where the magic happens)
  //NOTE: this is only one of several possible ways to code this
  uint64_t reassigned = threshold + 1;
  uint64_t localReassigned, specialReassigned = 0;
  const uint64_t REASSIGN_LIMIT = 1000000000; // M level recalc
  uint64_t reassignedSinceRecalculation = 0; // M level recalc
  int cycles = 0;
  while (reassigned > threshold) {
    it_start_time = MPI_Wtime();  // at start of iteration

    cycles++;
    localReassigned = 0;

    // reassignLocalDataPoints(options, localReassigned, data,
    //                         meanSize, meanValues, belongsToMean);
    reassignLocalDataPointsMlevel(options, localReassigned, data, meanSize, meanSizeLocal,
                                  meanValues, meanValuesLocal, belongsToMean); // M level recalc
    
    it_reassign_time = MPI_Wtime();  // after reassignment

    // now recalculate means (by communication) and check whether we are done
    // recalcMeans(options, data, meanSize, meanSizeBuff,
    //             meanValues, meanValuesBuff, belongsToMean); 
    MPI_Allreduce(&localReassigned, &reassigned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // M level recalc
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
                      meanValues, meanValuesBuff, meanValuesLocal, belongsToMean); // M level recalc
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

    // // communicate whether all processes have converged
    // MPI_Allreduce(&localReassigned, &reassigned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

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
  if (mpi_rank == 0) {
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " init time " << init_time - start_time << " seconds " << std::endl; 
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " reassign time "
              << cum_reassign_time
              << " | average  "
              << cum_reassign_time / cycles << " seconds" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " update time "
              << cum_update_time 
              << " | average "
              << cum_update_time / cycles << " seconds" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " mean recalculation time "
              << cum_it_end_time
              << " | average "
              << cum_it_end_time / cycles << " seconds" << std::endl;
  }

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
  delete[] meanSizeLocal; // M level recalc
  delete[] meanValuesLocal[0]; // M level recalc
  delete[] meanValuesLocal; // M level recalc
}


const uint64_t REASSIGN_LIMIT = 1000000000;  //used for testing, should be set so no rounding errors occur, doubles should have approx 15 digit precision


template <typename D, typename B>
void kmeansIncremental(struct Options &options, const std::string &variant,
                       D *data, B *belongsToMean)
{
  double start_time, end_time;

  // Time the individual parts of computation
  double it_start_time, it_reassign_time, it_update_time, it_end_time;
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
  
  //execute the algorithm (this is where the magic happens)
  //NOTE: this is only one of several possible ways to code this
  
  uint64_t reassigned = threshold + 1;
  uint64_t localReassigned, specialReassigned = 0;
  uint64_t reassignedSinceRecalculation = 0;
  int cycles = 0;
  while (reassigned > threshold) {
    it_start_time = MPI_Wtime();  // at start of iteration

    cycles++;
    localReassigned = 0;

    reassignLocalDataPoints(options, localReassigned, data,
                            meanSize, meanValues, belongsToMean);
    
    it_reassign_time = MPI_Wtime();  // after reassignment

    // communicate whether all processes have converged
    MPI_Allreduce(&localReassigned, &reassigned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //TODO: usefull to add a break if no reassignments are made?
    
    // now communicate the new values for the means between processes
    reassignedSinceRecalculation += reassigned;
    if (reassignedSinceRecalculation > REASSIGN_LIMIT || cycles < 2) {
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
    
    it_update_time = MPI_Wtime();  // after updating the means

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
  if (mpi_rank == 0) {
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " reassign time "
              << cum_reassign_time
              << " | average  "
              << cum_reassign_time / cycles << " seconds" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " update time "
              << cum_update_time 
              << " | average "
              << cum_update_time / cycles << " seconds" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " mean recalculation time "
              << cum_it_end_time
              << " | average "
              << cum_it_end_time / cycles << " seconds" << std::endl;
  }

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

// M level recalc
// This version can be more efficient if the local mean size is used from the start instead of added later.
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

  int mean;
  for (uint64_t x = 0; x < numLocalDataPoints; x++) {
    mean = getMean(belongsToMean, x);
    meanSizeLocal[mean]++;
    for (int d = 0; d < dataDim; d++) {
      meanValuesLocal[mean][d] += getDataPoint(data, x, d);
    }
  }

  for (int mean = 0; mean < numMeans; mean++) {
    for (int d = 0; d < dataDim; d++) {
      meanValuesLocal[mean][d] = meanValuesLocal[mean][d] / meanSizeLocal[mean];
    }
  }
}

// M level recalc
// Reassign but also keep track of mean size for local data points
template <typename D, typename B>
static inline void reassignLocalDataPointsMlevel(struct Options &options,
                                           uint64_t &localReassigned,
                                           D *data,
                                           uint64_t *meanSize,
                                           uint64_t * meanSizeLocal, // M level recalc
                                           double **meanValues,
                                           double **meanValuesLocal, // M level recalc
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
            // meanValues[oldMean][d] = (meanValues[oldMean][d] * meanSize[oldMean] - getDataPoint(data, x, d)) / (meanSize[oldMean] - 1);
            // meanValues[m][d] = (meanValues[m][d] * meanSize[m] + getDataPoint(data, x, d)) / (meanSize[m] + 1);

            meanValuesLocal[oldMean][d] = (meanValuesLocal[oldMean][d] * meanSizeLocal[oldMean] - getDataPoint(data, x, d)) / (meanSizeLocal[oldMean] - 1); // M level recalc
            meanValuesLocal[m][d] = (meanValuesLocal[m][d] * meanSizeLocal[m] + getDataPoint(data, x, d)) / (meanSizeLocal[m] + 1); // M level recalc
          }
          meanSize[oldMean] -= 1;
          meanSize[m] += 1;
          setMean(belongsToMean, x, m);
          meanSizeLocal[oldMean] -= 1; // M level recalc
          meanSizeLocal[m] += 1; // M level recalc
          oldDistance = newDistance; 
        }
      }
    }
  }

  // if (mpi_rank == 0) {
  //   std::cout << "; REASSIGNED: " << count << std::endl;
  // }
}

// reassign each of the local data points, updating the mean values as we go
// meanSize is kept strongly consistent. Each update is directly synchronized.
// Alternative version of reassignLocalDataPoints function.
template <typename D, typename B>
static inline void reassignSyncSize(struct Options &options,
                                    uint64_t &localReassigned,
                                    D *data,
                                    uint64_t *meanSize,
                                    double **meanValues,
                                    B *belongsToMean)
{
  const int numMeans(options.numMeans);
  const uint64_t numLocalDataPoints(options.numLocalDataPoints);
  const int dataDim(options.dataDim);

  uint64_t *meanSizeChange = new uint64_t[options.numMeans];

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

          /* Communicate change and update meanSize */
          // std::fill(meanSizeChange, meanSizeChange + options.numMeans, 0);
          // meanSizeChange[oldMean] -= 1;
          // meanSizeChange[m] += 1;
          // MPI_Allreduce(MPI_IN_PLACE, meanSizeChange, options.numMeans, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
          // for (int i = 0; i < options.numMeans; i++) {
          //   meanSize[i] += meanSizeChange[i];
          // }

          // V2: Recalculate meanSize through communication
          for (int i = 0; i < options.numMeans; i++) {
            meanSizeChange[i] = meanSize[i];
          }
          MPI_Allreduce(MPI_IN_PLACE, meanSizeChange, options.numMeans, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
          for (int i = 0; i < options.numMeans; i++) {
            meanSize[i] = meanSizeChange[i] / mpi_size;
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

/* New derivation variants */
template <typename D, typename B>
void kmeansRecalcNoDependencies(struct Options &options, const std::string &variant,
                  D *data, B *belongsToMean)
{
  double start_time, end_time;

  // Time the individual parts of computation
  double it_start_time, it_reassign_time, it_update_time, it_end_time;
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
    reassignSyncSize(options, localReassigned, data,
                            meanSize, meanValues, belongsToMean);
    
    it_reassign_time = MPI_Wtime();  // after reassignment

    // NEW: ONLY COMMUNICATE meanValuesBuff
    // now recalculate means (by communication) and check whether we are done
    recalcMeansValuesOnly(options, data, meanValues, meanValuesBuff, belongsToMean);
    
    it_update_time = MPI_Wtime();  // after updating the means

    // NEW: disabled (TODO: check if needed)
    // //check to avoid empty clusters, if so, reassign a point to fill them
    // ensureNoEmptyClusters(options, specialReassigned,
    //                       data, meanSize, meanValues, belongsToMean); 
    // broadcastMeans(options, meanValues, meanSize);

    // NEW: DIVIDE BY SIZE AND RANK COUNT
    divideMeansValuesOnly(options, meanValues, meanSize);
    
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
  if (mpi_rank == 0) {
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " reassign time "
              << cum_reassign_time
              << " | average  "
              << cum_reassign_time / cycles << " seconds" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " update time "
              << cum_update_time 
              << " | average "
              << cum_update_time / cycles << " seconds" << std::endl;
    std::cout << "TIME EXP " << options.currentRun << ": " << printPreamble()
              << " mean recalculation time "
              << cum_it_end_time
              << " | average "
              << cum_it_end_time / cycles << " seconds" << std::endl;
  }

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