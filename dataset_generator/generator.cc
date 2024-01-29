/**
 * MPI kmeans dataset generator. 
 * This generator constructs a dataset of a set size,
 * including initial means to be supplied to the algorithm.
 * 
 * Author: Dennis Buurman
 * Based on work of Anne Hommelberg 'generateData.cc'
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>

#include <sys/stat.h>

bool isdir(const char *path)
{
  struct stat statbuf;
  if (stat(path, &statbuf) != 0)
    return false;

  return S_ISDIR(statbuf.st_mode) != 0;
}

/**
 * Generates points in randomly placed clusters using normal distribution.
 * Points are generated in set size and duplicated for larger sized databases.
 * 
 * @param data_points: number of points generated
 * @param clusters: number of clusters generated
 * @param dimension: dimension of the data_points (and clusters) generated
 * @param output_dir:
 *      - generated data points are written to "data.txt" in this directory
 *      - the intended membership for each point is written to a file named "intended_membership.txt"
 *      - the generated cluster centres with their standard deviation and size are written to a file named "generated_cluster_centres.txt"
 *      - the generated initial means are saved in the file "initial_means.txt"
 * 
 * Notes: 
 *   - clusters may overlap
 *   - points are uniform randomly assigned to a cluster
 *   - generated cluster centres are uniformly generated in the interval [CENTRE_MIN,CENTRE_MAX] and standard deviations between 1/16 of and 1/8 of the width of this interval.
 */
int main(int argc, char *argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " [seed] [data_points] [clusters] [dimension] [output_dir]" << std::endl;
    }

    const char *outdir = argv[5];
    if (!isdir(outdir)) {
        std::cerr << "Error: " << outdir << " is not a directory or does not exist." << std::endl;
        return 1;
    }

    int seed = atoi(argv[1]);
    int power = atoi(argv[2]);
    uint64_t numDataPoints = 1 << power;
    int numClusters = atoi(argv[3]);
    int dataDim = atoi(argv[4]);

    std::cout << "Generating " << numDataPoints << " data points of dimension " << dataDim << " in " << numClusters << " clusters using random seed " << seed << ". Writing output to " << outdir << std::endl;
    
    double clusterCenters[numClusters][dataDim];
    double clusterStd[numClusters];
    uint64_t clusterSize[numClusters];

    std::default_random_engine generator(seed);

    // TODO

    return 0;
}