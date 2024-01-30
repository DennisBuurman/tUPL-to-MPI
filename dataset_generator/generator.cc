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

const double CENTRE_MIN = 0.0;
const double CENTRE_MAX = 10.0;

// Helper functions
double ** allocate2dArray(int rows, int columns) {
    double **array = new double *[rows];
    for (unsigned int i = 0; i < rows; i++) {
        array[i] = new double[columns];
    }
    return array;
}

bool isdir(const char *path)
{
  struct stat statbuf;
  if (stat(path, &statbuf) != 0)
    return false;

  return S_ISDIR(statbuf.st_mode) != 0;
}

/**
 * Dataset information struct
 * @param seed: dataset init seed
 * @param size: dataset size: 2^size
 * @param numClusters: amount of clusters
 * @param dataDim: amount of dimensions
 */
struct Dataset {
    int seed;
    int size;
    int numClusters;
    int dataDim;
    uint64_t numDataPoints;

    double **clusterCenters;  // 2d array of cluster centers
    double *clusterStd;  //
    uint64_t *clusterSize;

    Dataset(int seed_, int size_, int numClusters_, int dataDim_) {
        seed = seed_;
        size = size_;
        numClusters = numClusters_;
        dataDim = dataDim_;
        numDataPoints = 1 << size;

        // Allocate arrays
        clusterCenters = allocate2dArray(numClusters, dataDim);
        clusterStd = new double[numClusters];
        clusterSize = new uint64_t[numClusters];
    }

    ~Dataset() {
        // Deallocate arrays
        for (unsigned int i = 0; i < numClusters; i++) {
            delete[] clusterCenters[i];
        }
        delete[] clusterCenters;
        delete[] clusterStd;
        delete[] clusterSize;
    }
};

// Dataset functions
void generate_dataset(Dataset &data) {
    // TODO: random dataset of size 2^24
}

void upscale_dataset() {
    // TODO: upscale provided dataset from 2^x to 2^(x+y)
}

void downscale_dataset() {
    // TODO: downscale provided dataset from 2^x to 2^(x-y)
}

void generate_initial_means() {
    // TODO: generate (set) of initial means to be used during execution
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
    int size = atoi(argv[2]);
    int numClusters = atoi(argv[3]);
    int dataDim = atoi(argv[4]);

    Dataset *data = new Dataset(seed, size, numClusters, dataDim);

    std::cout << "Generating " << data->numDataPoints << " data points of dimension " << data->dataDim << " in " << data->numClusters << " clusters using random seed " << data->seed << ". Writing output to " << outdir << std::endl;
    
    // TODO

    return 0;
}