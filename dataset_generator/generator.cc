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
#include <vector>

#include <sys/stat.h>

const double CENTRE_MIN = 0.0;
const double CENTRE_MAX = 10.0;
const uint8_t DEFAULT_SIZE = 24;

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

    std::vector<std::vector<double>> clusterCenters;  // 2d array of cluster centers
    std::vector<double> clusterStd;  // array containing std dev of each cluster
    std::vector<uint64_t> clusterSize;  // array containing amount of points in each cluster
    std::vector<uint16_t> membership;  // array denoting intended cluster of datapoint at same index
    std::vector<std::vector<double>> datapoints;  // array containing the datapoints

    Dataset(const int seed_, const int size_, const int numClusters_, const int dataDim_) {
        seed = seed_;
        size = size_;
        numClusters = numClusters_;
        dataDim = dataDim_;
        numDataPoints = 1 << size;

        clusterSize.insert(clusterSize.begin(), numClusters, 0);
    }

    ~Dataset() {

    }
};

/**
 * Uniform randomly generate the means of the clusters in the interval [CENTRE_MIN, CENTRE_MAX]
 * @param dataset: Dataset object containing the dataset arrays and data
 * @param generator: reference to generator object used in sampling
*/
void generate_cluster_means(Dataset &dataset, std::default_random_engine &generator) {
    std::uniform_real_distribution<> centerdist(CENTRE_MIN,CENTRE_MAX);
    std::vector<double> center;

    for (unsigned int i = 0; i < dataset.numClusters; i++) {
        center.clear();
        for (unsigned int d = 0; d < dataset.dataDim; d++) {
            center.push_back(centerdist(generator));
        }
        dataset.clusterCenters.push_back(center);
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
        dataset.clusterStd.push_back(stddist(generator));
    }
}

/**
 * Create a dataset using the provided dataset object and random generator 
 * @param dataset: Dataset object containing the dataset arrays and data
 * @param generator: generator object used in sampling
*/
void generate_dataset(Dataset &dataset, std::default_random_engine &generator) {
    int cluster; // cluster number to put generated data point in
    std::vector<double> point; // vector containing point
    std::uniform_int_distribution<> clusterdist(0, dataset.numClusters-1);

    for (uint64_t i = 0; i < dataset.numDataPoints; i++) {
        // Randomly choose a cluster
        cluster = clusterdist(generator);
        dataset.membership.push_back(cluster);
        dataset.clusterSize[cluster]++;

        // Generate point
        point.clear();
        for (unsigned int d = 0; d < dataset.dataDim; d++) {
            std::normal_distribution<> dist(dataset.clusterCenters[cluster][d], dataset.clusterStd[cluster]);
            point.push_back(dist(generator));
        }
        dataset.datapoints.push_back(point);
    }
}

/**
 * Upscale provided dataset from 2^x to 2^(x+y)
 * @param dataset: Dataset object containing the dataset arrays and data
 * @param y: power increase of dataset size
*/
void upscale_dataset(Dataset &dataset, const int y) {
    // Allocate new datapoints
    for (int i = 0; i < y; i++) {
        // Duplicate dataset
        for (uint64_t j = 0; j < dataset.numDataPoints; j++) {
            dataset.clusterCenters.push_back(dataset.clusterCenters[j]);
            dataset.clusterStd.push_back(dataset.clusterStd[j]);
            dataset.membership.push_back(dataset.membership[j]);
            dataset.datapoints.push_back(dataset.datapoints[j]);
        }
        // Update dataset sizes
        dataset.size++;
        dataset.numDataPoints <<= 1;
    }
}

/**
 * Downscale provided dataset from 2^x to 2^(x-y)
 * @param dataset: Dataset object containing the dataset arrays and data
 * @param y: power decrease of dataset size
*/
void downscale_dataset(Dataset &dataset, const int y) {
    std::cout << "Not implemented yet!" << std::endl;
    // Question: how to reduce size below DEFAULT_SIZE?
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
    const int size = atoi(argv[2]);  // TODO: accept multiple sizes and generate in one go
    const int numClusters = atoi(argv[3]);
    const int dataDim = atoi(argv[4]);
    const int size_difference = size - DEFAULT_SIZE;

    std::default_random_engine generator(seed);
    Dataset *dataset = new Dataset(seed, DEFAULT_SIZE, numClusters, dataDim);

    std::cout << "Generating " << dataset->numDataPoints << " data points of dimension " << dataset->dataDim << " in " << dataset->numClusters << " clusters using random seed " << dataset->seed << ". Writing output to " << outdir << std::endl;
    
    // Generate cluster means and standard deviations
    generate_cluster_means(*dataset, generator);
    generate_cluster_std_devs(*dataset, generator);

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
    if (size_difference > 0) {
        upscale_dataset(*dataset, size_difference);
    } else if (size_difference < 0) {
        downscale_dataset(*dataset, size_difference);
    }

    // Write dataset to files
    write_data(*dataset, outdir);
    write_cluster_centers(*dataset, outdir);

    delete dataset;

    return 0;
}