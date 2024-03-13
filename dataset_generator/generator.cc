/**
 * MPI kmeans dataset generator source file. 
 * This generator constructs a dataset of a set size,
 * including initial means to be supplied to the algorithm.
 * 
 * Author: Dennis Buurman
 * Based on work of Anne Hommelberg 'generateData.cc'
*/

#include "generator.h"

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

Dataset::Dataset(const int seed_, const int size_, const int numClusters_, const int dataDim_, const char *outdir_) {
    seed = seed_;
    size = size_;
    numClusters = numClusters_;
    dataDim = dataDim_;
    numDataPoints = 1 << size;
    outdir = outdir_;
    generator.seed(seed);
    mean_generator.seed(seed);

    clusterSize.insert(clusterSize.begin(), numClusters, 0);
}

Dataset::~Dataset() {

}

void Dataset::print() const {
    unsigned int limit = 10;
    std::cout << "********************" << std::endl;
    std::cout << "Seed:      " << seed << std::endl;
    std::cout << "Size:      " << size << " (" << numDataPoints << ")" << std::endl;
    std::cout << "Clusters:  " << numClusters << std::endl;
    std::cout << "Dimension: " << dataDim << std::endl;
    std::cout << std::endl;
    std::cout << "Cluster Centres: " << "(" << clusterCenters.size() << ")" << std::endl;
    for (unsigned int i = 0; i < clusterCenters.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": (";
        for (unsigned int j = 0; j < clusterCenters[i].size(); j++) {
           std::cout << clusterCenters[i][j] << ", ";
        }
        std::cout << ")" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Cluster Standard Deviations: " << "(" << clusterStd.size() << ")" << std::endl;
    for (unsigned int i = 0; i < clusterStd.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": " << clusterStd[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Cluster Sizes: " << "(" << clusterSize.size() << ")" << std::endl;
    for (unsigned int i = 0; i < clusterSize.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": " << clusterSize[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Memberships: " << "(" << membership.size() << ")" << std::endl;
    for (unsigned int i = 0; i < membership.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": " << membership[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Datapoints: " << "(" << datapoints.size() << ")" << std::endl;
    for (unsigned int i = 0; i < datapoints.size(); i++) {
        if (i >= limit) 
            break;
        std::cout << i << ": (";
        for (unsigned int j = 0; j < datapoints[i].size(); j++) {
            std::cout << datapoints[i][j] << ", ";
        }
        std::cout << ")" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Initial means: " << "(" << initialMeans.size() << ")" << std::endl;
    for (unsigned int i = 0; i < initialMeans.size(); i++) {
        std::cout << "Mean set " << i << ": [";
        for (unsigned int j = 0; j < initialMeans[i].size(); j++) {
            std::cout << "(";
            for (unsigned int k = 0; k < datapoints[initialMeans[i][j]].size(); k++) {
                std::cout << datapoints[initialMeans[i][j]][k] << ", ";
            }
            std::cout << "), ";
        }
        std::cout << "]" << std::endl;
    }

    std::cout << "********************" << std::endl;
}

void Dataset::generate_cluster_means() {
    std::uniform_real_distribution<> centerdist(CENTRE_MIN,CENTRE_MAX);
    std::vector<double> center;

    for (unsigned int i = 0; i < numClusters; i++) {
        center.clear();
        for (unsigned int d = 0; d < dataDim; d++) {
            center.push_back(centerdist(generator));
        }
        clusterCenters.push_back(center);
    }
}

void Dataset::generate_cluster_std_devs() {
    std::uniform_real_distribution<> stddist((CENTRE_MAX-CENTRE_MIN)/16.0,(CENTRE_MAX-CENTRE_MIN)/8.0);

    for (unsigned int i = 0; i < numClusters; i++) {
        clusterStd.push_back(stddist(generator));
    }
}

void Dataset::generate_dataset() {
    int cluster; // cluster number to put generated data point in
    std::vector<double> point; // vector containing a point
    std::uniform_int_distribution<> clusterdist(0, numClusters-1);

    for (uint64_t i = 0; i < numDataPoints; i++) {
        // Randomly choose a cluster
        cluster = clusterdist(generator);
        membership.push_back(cluster);
        clusterSize[cluster]++;

        // Generate point
        point.clear();
        for (unsigned int d = 0; d < dataDim; d++) {
            std::normal_distribution<> dist(clusterCenters[cluster][d], clusterStd[cluster]);
            point.push_back(dist(generator));
        }
        datapoints.push_back(point);
    }
}

void Dataset::upscale_dataset(const unsigned int y) {
    std::cout << "Upscaling dataset from size " << size << " to size " << size + y << std::endl;
    // Allocate new datapoints
    for (unsigned int i = 0; i < y; i++) {
        membership.resize(2*numDataPoints);
        std::copy_n(membership.begin(), numDataPoints, membership.begin() + numDataPoints);
        datapoints.resize(2*numDataPoints);
        std::copy_n(datapoints.begin(), numDataPoints, datapoints.begin() + numDataPoints);
        std::transform(clusterSize.begin(), clusterSize.end(), clusterSize.begin(), std::bind1st(std::multiplies<uint64_t>(), 2));
        // Update dataset size
        size++;
        numDataPoints <<= 1;
    }
    if (datapoints.size() != numDataPoints) {
        std::cerr << "ERROR: amount of datapoints " << datapoints.size() << " differs from expected amount " << numDataPoints << std::endl;
    }
}

double Dataset::distance(const std::vector<double> x, const std::vector<double> y) const {
    double d = 0;
    if (x.size() != y.size()) {
        std::cerr << "WARNING: points x (" << x.size() << ") and y (" << y.size() << ") are of different dimension. " << std::endl;
        return 0;
    }
    for (unsigned int i = 0; i < x.size(); i++) {
        d += (x[i] - y[i]) * (x[i] - y[i]);
    }
    // return sqrt(d);
    return d;
}

void Dataset::generate_initial_means() {
    uint64_t centroid_index;
    double d, d_min = DOUBLE_MAX, d_max = -1; // distance variables
    std::vector<uint64_t> means; // vector for a set of initial means
    std::vector<double> point, centroid;

    // Randomly select first centroid (uniformly)
    // std::random_device rand_dev;
    // std::mt19937 g(rand_dev());
    std::uniform_int_distribution<uint64_t> dist(0, numDataPoints-1);
    centroid_index = dist(mean_generator);
    means.push_back(centroid_index);

    // Select remaining k - 1 initial cluster centroids
    for (unsigned int k = 0; k < numClusters - 1; k++) {
        d_max = -1;
        
        // Select point with maximum distance from previous centroids as next centroid
        for (uint64_t i = 0; i < datapoints.size(); i++) {
            point = datapoints[i];
            d_min = DOUBLE_MAX;
            if (!in_vector(i, means)) {
                // Calculate minimum distance between point i and previously chosen centroids
                for (unsigned int j = 0; j < means.size(); j++) {
                    centroid = datapoints[means[j]];
                    d = distance(point, centroid);
                    d = std::min(d, d_min);
                }
                // Select max distance from closest centroid distance for each point
                if (d_min > d_max) {
                    centroid_index = i;
                    d_max = d_min;
                }
            }
        }

        // Add new centroid to initial means
        means.push_back(centroid_index);
    }
    initialMeans.push_back(means);
}

int Dataset::write_data() const {
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

    for (uint64_t i = 0; i < numDataPoints; i++) {
        // Write point i
        datafile << i + 1 << " ";
        for (unsigned int d = 0; d < dataDim; d++) {
            datafile << datapoints[i][d] << " ";
        }
        datafile << std::endl;
        // Write intended membership of point i
        memberfile << i + 1 << " " << membership[i] << std::endl;
    }

    datafile.close();
    memberfile.close();
    return 0;
}

int Dataset::write_cluster_centers() const {
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
    for (unsigned int i = 0; i < numClusters; i++) {
        centerfile << i << " ";
        for (unsigned int d = 0; d < dataDim; d++) {
            centerfile << clusterCenters[i][d] << " ";
        }
        centerfile << ", size = " << clusterSize[i] << " , std = " << clusterStd[i] << std::endl;
    }

    centerfile.close();
    return 0;
}

int Dataset::write_initial_means() const {
    std::stringstream ss_means;
    std::ofstream meansfile;
    std::vector<double> point;

    for (int i = 0; i < MEAN_SETS; i++) {
        // Open initial means file corresponding to index + 1
        ss_means.clear();
        ss_means.str(std::string());
        ss_means << outdir << "/initial_means_" << i+1 << ".txt";
        meansfile.open(ss_means.str());

        // Check if file is open
        if (!meansfile.is_open()) {
            std::cerr << "Error opening initial means file " << i << " in " << outdir << ", terminating..." << std::endl;
            return 2;
        }

        // Write initial mean set to file
        for (unsigned int j = 0; j < initialMeans[i].size(); j++) {
            meansfile << j << " ";
            point = datapoints[initialMeans[i][j]];
            for (unsigned int k = 0; k < point.size(); k++) {
                meansfile << point[k] << " ";
            }
            meansfile << std::endl;
        }
        
        // Close file and reset flags
        meansfile.close();
        meansfile.clear();
    }

    return 0;
}

int Dataset::write_metadata() const {
    // Open cluster center file
    std::stringstream ss_metadata;
    ss_metadata << outdir << "/metadata.md";
    std::ofstream metadatafile(ss_metadata.str());

    // Check if file is open
    if (!metadatafile.is_open()) {
        std::cerr << "Error opening metadata file in " << outdir << ", terminating..." << std::endl;
        return 2;
    }

    // Write metadata from dataset and program settings
    metadatafile << "# Metadata file for k-means dataset generator" << std::endl;
    metadatafile << "seed: " << seed << std::endl;
    metadatafile << "size: " << size << std::endl;
    metadatafile << "datapoints: " << numDataPoints << std::endl;
    metadatafile << "dimension: " << dataDim << std::endl;
    metadatafile << "clusters: " << numClusters << std::endl;
    metadatafile << "default_size: " << MINIMUM_SIZE << std::endl;
    metadatafile << "mean_sets: " << MEAN_SETS << std::endl;

    metadatafile.close();
    return 0;
}

void Dataset::write() const {
    write_data();
    write_cluster_centers();
    write_initial_means();
    write_metadata();
}
