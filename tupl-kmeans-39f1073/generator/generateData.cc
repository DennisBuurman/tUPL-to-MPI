/*
 * MPI K-means data generator.
 *
 * Author: Anne Hommelberg
 */


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <random>
#include <algorithm>

#include <sys/stat.h>


using namespace std;

const double CENTRE_MIN = 0.0;
const double CENTRE_MAX = 10.0;



bool isdir(const char *path)
{
  struct stat statbuf;
  if (stat(path, &statbuf) != 0)
    return false;

  return S_ISDIR(statbuf.st_mode) != 0;
}

/**
 * Generates points in randomly placed clusters using normal distribution
 * Usage: exec_name [data points] [clusters] [dimension] [output dir]
 *
 *  @param data points: number of points generated
 *  @param clusters: number of clusters generated
 *  @param dimension: dimension of the data points (and clusters) generated
 *  @param output dir: 
 *          - generated data points are written to "data.txt" in this directory
 *          - the intended membership for each point is written to a file named "intended_membership.txt"
 *          - the generated cluster centres with their standard deviation and size are written to a file named "generated_cluster_centres.txt"
 *
 * Notes: 
 *   - clusters may overlap
 *   - points are uniform randomly assigned to a cluster
 *   - generated cluster centres are uniformly generated in the interval [CENTRE_MIN,CENTRE_MAX] and standard deviations between 1/16 of and 1/8 of the width of this interval.
 */
int main (int argc, char * argv[]) {
  if (argc != 6) {
    cerr << "Usage: " << argv[0] << " [seed] [data points] [clusters] [dimension] [output dir]" << endl;
    return 1;
  }

  const char *outdir = argv[5];
  if (!isdir(outdir)) {
    cerr << "Error: " << outdir << " is not a directory or does not exist." << endl;
    return 1;
  }

  int seed = atoi(argv[1]);
  int power = atoi(argv[2]);
  uint64_t numDataPoints = 2;
  for (int i = 1; i < power; i++) {
    numDataPoints = numDataPoints * 2;
  }
  int numClusters = atoi(argv[3]);
  int dataDim = atoi(argv[4]);
  
  cout << "Generating " << numDataPoints << " data points of dimension " << dataDim << " in " << numClusters << " clusters using random seed " << seed << ". Writing output to " << outdir << endl;
  
  double clusterCenters[numClusters][dataDim];
  double clusterStd[numClusters];
  uint64_t clusterSize[numClusters];
  
  default_random_engine generator(seed);
  
  //uniform randomly generate the means of the clusters in the interval [CENTRE_MIN, CENTRE_MAX]
  uniform_real_distribution<> centerdist(CENTRE_MIN,CENTRE_MAX);
  for (int i = 0; i < numClusters; i++) {
    for (int d = 0; d < dataDim; d++) {
      clusterCenters[i][d] = centerdist(generator);
    }
  }
  
  //uniform randomly generate the standard deviation of the clusters between 1/16 and 1/8 of the width of the interval [CENTRE_MIN, CENTRE_MAX]
  uniform_real_distribution<> stddist((CENTRE_MAX-CENTRE_MIN)/16.0,(CENTRE_MAX-CENTRE_MIN)/8.0);
  for (int i = 0; i < numClusters; i++) {
    clusterStd[i] = stddist(generator);
  }
  
  //initialize the sizes on zero
  fill(clusterSize, clusterSize + numClusters, 0);
  
  //some output for the user
  cout << endl << "Generating " << numDataPoints << " clustered data points, with clusters:" << endl;
  for (int i = 0; i < numClusters; i++) {
    cout << "  " << i << ": (";
    for (int d = 0; d < dataDim; d++) {
      cout << clusterCenters[i][d];
      if (d != dataDim-1) cout << ",";
    }
    cout << "), std: " << clusterStd[i] << endl;
  }
  cout << endl;
  
  //open the output files for the data points and their intended membership
  stringstream ss_data;
  ss_data << outdir << "/data.txt";
  ofstream datafile(ss_data.str());

  stringstream ss_members;
  ss_members << outdir << "/intended_membership.txt";
  ofstream memberfile(ss_members.str());
  if (!datafile.is_open() || !memberfile.is_open()) {
    cerr << "Error opening output files in " << outdir << ", terminating..." << endl;
    return 2;
  }
  
  //generate the data and write it to the files as we go
  uniform_int_distribution<> clusterdist(0,numClusters-1);
  for (uint64_t i = 0; i < numDataPoints; i++) {
    //randomly choose a cluster
    int cluster = clusterdist(generator);
    memberfile << i + 1 << " " << cluster << endl;
    clusterSize[cluster]++;
    
    //now use normal distribution with specified mean and std to generate point
    datafile << i + 1 << " ";
    for (int d = 0; d < dataDim; d++) {
      normal_distribution<> dist(clusterCenters[cluster][d],clusterStd[cluster]);
      datafile << dist(generator) << " ";
    }
    datafile << endl;
  }
  
  datafile.close();
  memberfile.close();
  
  //print the generated clusters to file
  stringstream ss_centres;
  ss_centres << argv[5] << "/generated_cluster_centres.txt";
  ofstream centerfile(ss_centres.str());

  for (int i = 0; i < numClusters; i++) {
    centerfile << i << " ";
    for (int d = 0; d < dataDim; d++) {
      centerfile << clusterCenters[i][d] << " ";
    }
    centerfile << ", size = " << clusterSize[i] << " , std = " << clusterStd[i] << endl;
  }
  centerfile.close();
  
  return 0;
}
