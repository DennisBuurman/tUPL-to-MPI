/**
 * MPI K-means implementation 
 * Author: Dennis Buurman
 * Code is based on work of Anne Hommelberg (MPI_Kmeans.cpp) 
 */

#include "algorithm.h"

#include <iostream>


// Executes k-means algorithm 
// NOTE: the datapoints are distributed across processes, process 0 is required to receive at least k data points, k should be selected accordingly, this is not checked!
int kmeans(struct Options &options) {
  
  //read input from file
  double ** data = readDataAsDoubles(options);
  int * belongsToMean = new int[options.numLocalDataPoints];

  for (options.currentRun = 1; options.currentRun <= options.numRuns;
       options.currentRun++)
    kmeansRecalcNoUpdates(options, "own_im", data, belongsToMean);
  
  delete[] data[0];
  delete[] data;
  delete[] belongsToMean;
  return 0;
}

int main (int argc, char ** argv) {
  return runVariant(argc, argv, kmeans);
}