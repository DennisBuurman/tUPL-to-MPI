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
#include <limits>

#include <sys/stat.h>

const double CENTRE_MIN = 0.0;
const double CENTRE_MAX = 10.0;
const double DOUBLE_MAX = std::numeric_limits<double>::max();

const uint8_t MINIMUM_SIZE = 24;
const uint8_t MEAN_SETS = 10;

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

    std::vector<std::vector<double>> clusterCenters;  // 2d vector of cluster centers
    std::vector<double> clusterStd;  // vector containing std dev of each cluster
    std::vector<uint64_t> clusterSize;  // vector containing amount of points in each cluster
    std::vector<uint32_t> membership;  // vector denoting intended cluster of datapoint at same index
    std::vector<std::vector<double>> datapoints;  // vector containing the datapoints
    std::vector<std::vector<uint64_t>> initialMeans;  // vector containing indexes of initial means

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
 * Prints a given dataset
 * @param data: dataset to print
*/
void print_dataset(Dataset &data) {
    unsigned int limit = 10;
    std::cout << "********************" << std::endl;
    std::cout << "Seed:      " << data.seed << std::endl;
    std::cout << "Size:      " << data.size << " (" << data.numDataPoints << ")" << std::endl;
    std::cout << "Clusters:  " << data.numClusters << std::endl;
    std::cout << "Dimension: " << data.dataDim << std::endl;
    std::cout << std::endl;
    std::cout << "Cluster Centres: " << "(" << data.clusterCenters.size() << ")" << std::endl;
    for (unsigned int i = 0; i < data.clusterCenters.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": (";
        for (unsigned int j = 0; j < data.clusterCenters[i].size(); j++) {
           std::cout << data.clusterCenters[i][j] << ", ";
        }
        std::cout << ")" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Cluster Standard Deviations: " << "(" << data.clusterStd.size() << ")" << std::endl;
    for (unsigned int i = 0; i < data.clusterStd.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": " << data.clusterStd[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Cluster Sizes: " << "(" << data.clusterSize.size() << ")" << std::endl;
    for (unsigned int i = 0; i < data.clusterSize.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": " << data.clusterSize[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Memberships: " << "(" << data.membership.size() << ")" << std::endl;
    for (unsigned int i = 0; i < data.membership.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": " << data.membership[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Datapoints: " << "(" << data.datapoints.size() << ")" << std::endl;
    for (unsigned int i = 0; i < data.datapoints.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": (";
        for (unsigned int j = 0; j < data.datapoints[i].size(); j++) {
            std::cout << data.datapoints[i][j] << ", ";
        }
        std::cout << ")" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Initial means: " << "(" << data.initialMeans.size() << ")" << std::endl;
    for (unsigned int i = 0; i < data.initialMeans.size(); i++) {
        std::cout << "Mean set " << i << ": [";
        for (unsigned int j = 0; j < data.initialMeans[i].size(); j++) {
            std::cout << "(";
            for (unsigned int k = 0; k < data.datapoints[data.initialMeans[i][j]].size(); k++) {
                std::cout << data.datapoints[data.initialMeans[i][j]][k] << ", ";
            }
            std::cout << "), ";
        }
        std::cout << "]" << std::endl;
    }

    std::cout << "********************" << std::endl;
}

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
    std::vector<double> point; // vector containing a point
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
    std::cout << "Upscaling dataset from size " << dataset.size << " to size " << dataset.size + y << std::endl;
    // Allocate new datapoints
    for (int i = 0; i < y; i++) {
        // Duplicate dataset
        for (uint64_t j = 0; j < dataset.numDataPoints; j++) {
            dataset.clusterCenters.push_back(dataset.clusterCenters[j]);
            dataset.clusterStd.push_back(dataset.clusterStd[j]);
            dataset.membership.push_back(dataset.membership[j]);
            dataset.datapoints.push_back(dataset.datapoints[j]);
            dataset.clusterSize[dataset.membership[j]]++;
        }
        // Update dataset sizes
        dataset.size++;
        dataset.numDataPoints <<= 1;
    }
}

/**
 * Euclidean distance between two points.
 * Dimensions between points x and y need to be identical.
 * @param x: point 1
 * @param y: point 2
*/
double distance(std::vector<double> x, std::vector<double> y) {
    double d = 0;
    if (x.size() != y.size()) {
        std::cerr << "WARNING: points x (" << x.size() << ") and y (" << y.size() << ") are of different dimension. " << std::endl;
        return 0;
    }
    for (unsigned int i = 0; i < x.size(); i++) {
        d += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(d);
}

/**
 * Generate (set) of initial means to be used during execution.
 * This function uses the k-means++ method to generate initial means.
 * In short, iteratively selects the point with maximum distance from previously selected centroids.
 * First centroid is randomly selected.
 * Selected centroids are used as initial means.
 * Initial means are saved in initial_means vector of dataset.
 * Based on implementation from https://www.geeksforgeeks.org/ml-k-means-algorithm/.
 * 
 * @param dataset: Dataset object containing the dataset arrays and data
*/
void generate_initial_means(Dataset &dataset) {
    uint64_t centroid_index;
    double d, d_min = DOUBLE_MAX, d_max = -1; // distance between 2 points
    std::vector<uint64_t> means;

    // Randomly select first centroid
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<uint64_t> dist(0, dataset.numDataPoints);
    centroid_index = dist(generator);
    means.push_back(centroid_index);

    // Select remaining k - 1 initial cluster centroids
    for (unsigned int k = 0; k < dataset.numClusters - 1; k++) {
        // Select point with maximum distance from previous centroids as next centroid
        for (unsigned int i = 0; i < dataset.datapoints.size(); i++) { 
            // Calculate minimum distance from nearest, previously chosen centroids
            d_min = DOUBLE_MAX;
            for (unsigned int j = 0; j < dataset.initialMeans.size(); j++) {
                d = distance(dataset.datapoints[i], dataset.datapoints[means[j]]);
                if (d < d_min) {
                    d_min = d;
                }
            }
            // Select max distance from closest centroid distance for each point
            if (d_min > d_max) {
                centroid_index = i;
                d_max = d_min;
            }
        }
        // Add new centroid to initial means
        means.push_back(centroid_index);
    }
    dataset.initialMeans.push_back(means);
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
 * Write generated means to files
 * @param dataset: dataset object containing data to write
 * @param outdir: output file directory
*/
int write_initial_means(Dataset &dataset, const char *outdir) {
    std::stringstream ss_means;
    std::ofstream meansfile;

    for (int i = 0; i < MEAN_SETS; i++) {
        // Open initial means file corresponding to index + 1
        ss_means.clear();
        ss_means << outdir << "/initial_means" << i+1 << ".txt";
        meansfile.open(ss_means.str());

        // Check if file is open
        if (!meansfile.is_open()) {
            std::cerr << "Error opening cluster centre file in " << outdir << ", terminating..." << std::endl;
            return 2;
        }

        // Write cluster centres to file
        // TODO
        
        meansfile.close();
        meansfile.clear();
    }

    return 0;
}

/**
 * Function for testing the dataset generations.
 * Does running this function 2 times result in identical datasets? (should be yes)
*/
void test() {
    int seed = 971, size = 8, clusters = 4, dim = 4, upscale = 1, mean_sets = 3;

    std::default_random_engine g(seed);
    Dataset *d = new Dataset(seed, size, clusters, dim);
    generate_cluster_means(*d, g);
    generate_cluster_std_devs(*d, g);
    generate_dataset(*d, g);
    upscale_dataset(*d, upscale);
    for (int i = 0; i < mean_sets; i++) {
        generate_initial_means(*d);
    }

    print_dataset(*d);
    delete d;
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
    const int size_difference = size - MINIMUM_SIZE;

    // // DEBUGGING ONLY
    // test();
    // return 0;

    std::default_random_engine generator(seed);
    Dataset *dataset = new Dataset(seed, MINIMUM_SIZE, numClusters, dataDim);

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
        std::cerr << "ERROR: requested size " << size << " smaller than minimum size " << MINIMUM_SIZE << std::endl;
    }

    // Generate sets of initial means that can be supplied to algorithms
    for (int i = 0; i < MEAN_SETS; i++) {
        generate_initial_means(*dataset);
    }

    if (size != dataset->size) {
        std::cerr << "ERROR: size " << size << " smaller than dataset size " << dataset->size << std::endl;
        return 1;
    }

    // Write dataset to files
    write_data(*dataset, outdir);
    write_cluster_centers(*dataset, outdir);
    write_initial_means(*dataset, outdir);

    delete dataset;

    return 0;
}