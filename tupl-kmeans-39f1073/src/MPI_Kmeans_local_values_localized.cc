/**
 * MPI K-means implementation 
 * Author: Dennis Buurman
 * Code is based on work of Anne Hommelberg (MPI_Kmeans_localized.cpp) 
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
    kmeansRecalcMlevel(options, "own_m_loc", data, data);
  
  delete[] data;
  return 0;
}

int main (int argc, char ** argv) {
  return runVariant(argc, argv, kmeans);
}
