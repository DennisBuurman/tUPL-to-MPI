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
const uint8_t DEFAULT_CLUSTER_SIZE = 24;

double ** allocate2dArray(const unsigned int rows, const unsigned int columns) {
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
    unsigned int numClusters;
    unsigned int dataDim;
    uint64_t numDataPoints;

    double **clusterCenters;  // 2d array of cluster centers
    double *clusterStd;  // array containing std dev of each cluster
    uint64_t *clusterSize;  // array containing amount of points in each cluster
    uint16_t *membership;  // array denoting intended cluster of datapoint at same index
    double **datapoints;  // array containing the datapoints

    Dataset(const int seed_, const int size_, const int numClusters_, const int dataDim_) {
        seed = seed_;
        size = size_;
        numClusters = numClusters_;
        dataDim = dataDim_;
        numDataPoints = 1 << size;

        // Allocate arrays
        clusterCenters = allocate2dArray(numClusters, dataDim);
        clusterStd = new double[numClusters];
        clusterSize = new uint64_t[numClusters];
        membership = new uint16_t[numDataPoints];
        datapoints = allocate2dArray(numDataPoints, dataDim);
    }

    ~Dataset() {
        // Deallocate arrays
        for (unsigned int i = 0; i < numClusters; i++) {
            delete[] clusterCenters[i];
        }
        delete[] clusterCenters;
        delete[] clusterStd;
        delete[] clusterSize;
        delete[] membership;
        for (unsigned int i = 0; i < numDataPoints; i++) {
            delete[] datapoints[i];
        }
    }
};

/**
 * Uniform randomly generate the means of the clusters in the interval [CENTRE_MIN, CENTRE_MAX]
 * @param dataset: Dataset object containing the dataset arrays and data
 * @param generator: reference to generator object used in sampling
*/
void generate_cluster_means(Dataset &dataset, std::default_random_engine &generator) {
    std::uniform_real_distribution<> centerdist(CENTRE_MIN,CENTRE_MAX);
    for (unsigned int i = 0; i < dataset.numClusters; i++) {
        for (unsigned int d = 0; d < dataset.dataDim; d++) {
            dataset.clusterCenters[i][d] = centerdist(generator);
        }
    }
}

/**
 * Uniform randomly generate the standard deviation of the clusters between 1/16 and 1/8 of the width of the interval [CENTRE_MIN, CENTRE_MAX]
 * @param dataset: Dataset object containing the dataset arrays and data
 * @param generator: generator object used in sampling
*/
void generate_cluster_std_devs(Dataset &dataset, std::default_random_engine &generator) {
    std::uniform_real_distribution<> stddist((CENTRE_MAX-CENTRE_MIN)/16.0,(CENTRE_MAX-CENTRE_MIN)/8.0);
    for (unsigned int i = 0; i < dataset.numClusters; i++) {
        dataset.clusterStd[i] = stddist(generator);
    }
}

/**
 * Create a dataset using the provided dataset object and random generator 
 * @param dataset: Dataset object containing the dataset arrays and data
 * @param generator: generator object used in sampling
*/
void generate_dataset(Dataset &dataset, std::default_random_engine &generator) {
    int cluster; // cluster number to put generated data point in

    std::uniform_int_distribution<> clusterdist(0, dataset.numClusters-1);
    for (uint64_t i = 0; i < dataset.numDataPoints; i++) {
        // Randomly choose a cluster
        cluster = clusterdist(generator);
        dataset.membership[i] = cluster;
        dataset.clusterSize[cluster]++;

        // Generate point
        for (unsigned int d = 0; d < dataset.dataDim; d++) {
            std::normal_distribution<> dist(dataset.clusterCenters[cluster][d], dataset.clusterStd[cluster]);
            dataset.datapoints[i][d] = dist(generator);
        }
    }
}

/**
 * Upscale provided dataset from 2^x to 2^(x+y)
 * @param TODO
*/
void upscale_dataset(Dataset &dataset, const int y) {
    const int x = dataset.size;
    const uint64_t new_numDataPoints = 1 << x+y; // equivalent to 2^(x+y)
    const uint64_t difference = new_numDataPoints - dataset.numDataPoints;

    // Allocate new datapoints
    for (uint64_t i = dataset.numDataPoints; i < new_numDataPoints; i++) {

    }
}

/**
 * Downscale provided dataset from 2^x to 2^(x-y)
 * @param TODO
*/
void downscale_dataset(Dataset &dataset, const int y) {
    // TODO
}

/**
 * Generate (set) of initial means to be used during execution
 * @param TODO
*/
void generate_initial_means(Dataset &dataset) {
    // TODO
}

/**
 * Write data points and their intended membership to corresponding files
 * @param dataset: dataset object containing data to write
 * @param outdir: output file directory
*/
int write_data(Dataset &dataset, const char *outdir) {
    // Open data file
    std::stringstream ss_data;
    ss_data << outdir << "/data.txt";
    std::ofstream datafile(ss_data.str());

    // Open membership file
    std::stringstream ss_members;
    ss_members << outdir << "/intended_membership.txt";
    std::ofstream memberfile(ss_members.str());

    // Check if files are open
    if (!datafile.is_open() || !memberfile.is_open()) {
        std::cerr << "Error opening output files in " << outdir << ", terminating..." << std::endl;
        return 2;
    }

    for (uint64_t i = 0; i < dataset.numDataPoints; i++) {
        // Write point i
        datafile << i + 1 << " ";
        for (unsigned int d = 0; d < dataset.dataDim; d++) {
            datafile << dataset.datapoints[i][d] << " ";
        }
        datafile << std::endl;
        // Write intended membership of point i
        memberfile << i + 1 << " " << dataset.membership[i] << std::endl;
    }

    datafile.close();
    memberfile.close();
    return 0;
}

/**
 * Write generated cluster to file
 * @param dataset: dataset object containing data to write
 * @param outdir: output file directory
*/
int write_cluster_centers(Dataset &dataset, const char *outdir) {
    // Open cluster center file
    std::stringstream ss_centres;
    ss_centres << outdir << "/generated_cluster_centres.txt";
    std::ofstream centerfile(ss_centres.str());

    // Check if file is open
    if (!centerfile.is_open()) {
        std::cerr << "Error opening cluster centre file in " << outdir << ", terminating..." << std::endl;
        return 2;
    }

    // Write cluster centres to file
    for (unsigned int i = 0; i < dataset.numClusters; i++) {
        centerfile << i << " ";
        for (unsigned int d = 0; d < dataset.dataDim; d++) {
            centerfile << dataset.clusterCenters[i][d] << " ";
        }
        centerfile << ", size = " << dataset.clusterSize[i] << " , std = " << dataset.clusterStd[i] << std::endl;
    }

    centerfile.close();
    return 0;
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

    const int seed = atoi(argv[1]);
    const int size = atoi(argv[2]);
    const int numClusters = atoi(argv[3]);
    const int dataDim = atoi(argv[4]);

    std::default_random_engine generator(seed);
    Dataset *dataset = new Dataset(seed, DEFAULT_CLUSTER_SIZE, numClusters, dataDim);

    std::cout << "Generating " << dataset->numDataPoints << " data points of dimension " << dataset->dataDim << " in " << dataset->numClusters << " clusters using random seed " << dataset->seed << ". Writing output to " << outdir << std::endl;
    
    // Generate cluster means and standard deviations
    generate_cluster_means(*dataset, generator);
    generate_cluster_std_devs(*dataset, generator);
    std::fill(dataset->clusterSize, dataset->clusterSize + dataset->numClusters, 0); //initialize the sizes on zero

    // Some output for the user
    std::cout << std::endl << "Generating " << dataset->numDataPoints << " clustered data points, with clusters:" << std::endl;
    for (int i = 0; i < numClusters; i++) {
        std::cout << "  " << i << ": (";
        for (int d = 0; d < dataDim; d++) {
            std::cout << dataset->clusterCenters[i][d];
            if (d != dataDim-1) 
                std::cout << ",";
        }
        std::cout << "), std: " << dataset->clusterStd[i] << std::endl;
    }
    std::cout << std::endl;

    // Generate default dataset and upscale/downscale it
    generate_dataset(*dataset, generator);
    

    // Write dataset to files
    write_data(*dataset, outdir);
    write_cluster_centers(*dataset, outdir);

    delete dataset;

    return 0;
}