/*
 * MPI K-means implementation
 *
 * Author: Anne Hommelberg
 *
 * Code refactoring done by Kristian Rietveld, Leiden University.
 */

#include "algorithm.h"

#include <iostream>


// Executes k-means algorithm 
// NOTE: the datapoints are distributed across processes, process 0 is required to receive at least k data points, k should be selected accordingly, this is not checked!
int kmeans(struct Options &options) {
  
  defineMPIDataPoint();
  
  //read input from file
  DataPoint * data = readDataAsDataPoints(options);
  
  for (options.currentRun = 1; options.currentRun <= options.numRuns;
       options.currentRun++)
    kmeansRecalc(options, "own_loc", data, data);
  
  delete[] data;
  return 0;
}

int main (int argc, char ** argv) {
  return runVariant(argc, argv, kmeans);
}
