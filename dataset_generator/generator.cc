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

bool isdir(const char *path)
{
  struct stat statbuf;
  if (stat(path, &statbuf) != 0)
    return false;

  return S_ISDIR(statbuf.st_mode) != 0;
}

void generate_dataset() {
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
    int power = atoi(argv[2]);
    uint64_t numDataPoints = 1 << power;
    int numClusters = atoi(argv[3]);
    int dataDim = atoi(argv[4]);

    std::cout << "Generating " << numDataPoints << " data points of dimension " << dataDim << " in " << numClusters << " clusters using random seed " << seed << ". Writing output to " << outdir << std::endl;
    
    double clusterCenters[numClusters][dataDim];
    double clusterStd[numClusters];
    uint64_t clusterSize[numClusters];

    std::default_random_engine generator(seed);

    //uniform randomly generate the means of the clusters in the interval [CENTRE_MIN, CENTRE_MAX]
    std::uniform_real_distribution<> centerdist(CENTRE_MIN,CENTRE_MAX);
    for (int i = 0; i < numClusters; i++) {
        for (int d = 0; d < dataDim; d++) {
            clusterCenters[i][d] = centerdist(generator);
        }
    }

    //initialize the sizes on zero
    std::fill(clusterSize, clusterSize + numClusters, 0);

    //some output for the user
    std::cout << std::endl << "Generating " << numDataPoints << " clustered data points, with clusters:" << std::endl;
    for (int i = 0; i < numClusters; i++) {
        std::cout << "  " << i << ": (";
        for (int d = 0; d < dataDim; d++) {
            std::cout << clusterCenters[i][d];
            if (d != dataDim-1) 
                std::cout << ",";
        }
        std::cout << "), std: " << clusterStd[i] << std::endl;
    }
    std::cout << std::endl;

    //open the output files for the data points and their intended membership
    std::stringstream ss_data;
    ss_data << outdir << "/data.txt";
    std::ofstream datafile(ss_data.str());

    std::stringstream ss_members;
    ss_members << outdir << "/intended_membership.txt";
    std::ofstream memberfile(ss_members.str());
    if (!datafile.is_open() || !memberfile.is_open()) {
        std::cerr << "Error opening output files in " << outdir << ", terminating..." << std::endl;
        return 2;
    }

    //generate the data and write it to the files as we go
    std::uniform_int_distribution<> clusterdist(0,numClusters-1);
    for (uint64_t i = 0; i < numDataPoints; i++) {
        //randomly choose a cluster
        int cluster = clusterdist(generator);
        memberfile << i + 1 << " " << cluster << std::endl;
        clusterSize[cluster]++;
        
        //now use normal distribution with specified mean and std to generate point
        datafile << i + 1 << " ";
        for (int d = 0; d < dataDim; d++) {
            std::normal_distribution<> dist(clusterCenters[cluster][d],clusterStd[cluster]);
            datafile << dist(generator) << " ";
        }
        datafile << std::endl;
    }
    
    datafile.close();
    memberfile.close();

    //print the generated clusters to file
    std::stringstream ss_centres;
    ss_centres << argv[5] << "/generated_cluster_centres.txt";
    std::ofstream centerfile(ss_centres.str());

    for (int i = 0; i < numClusters; i++) {
        centerfile << i << " ";
        for (int d = 0; d < dataDim; d++) {
        centerfile << clusterCenters[i][d] << " ";
        }
        centerfile << ", size = " << clusterSize[i] << " , std = " << clusterStd[i] << std::endl;
    }
    centerfile.close();

    return 0;
}