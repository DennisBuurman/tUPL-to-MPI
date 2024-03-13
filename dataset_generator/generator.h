/**
 * MPI kmeans dataset generator header file. 
 * This generator constructs a dataset of a set size,
 * including initial means to be supplied to the algorithm.
 * 
 * Author: Dennis Buurman
 * Based on work of Anne Hommelberg 'generateData.cc'
*/

#ifndef __GENERATOR_H__
#define __GENERATOR_H__

#include <iostream>
#include <random>
#include <vector>
#include <limits>
#include <string>

const double CENTRE_MIN = 0.0;
const double CENTRE_MAX = 10.0;
const double DOUBLE_MAX = std::numeric_limits<double>::max();

const int MINIMUM_SIZE = 24; // minimum required size of the dataset
const int MEAN_SETS = 10; // amount of mean sets to generate

/**
 * Check if uint64_t x is in vector.
 * @param x: element to check
 * @param v: vector
*/
inline bool in_vector(uint64_t x, std::vector<uint64_t> v) {
    for (uint64_t i = 0; i < v.size(); i++) {
        if (x == v[i]) {
            return true;
        }
    }
    return false;
}

/**
 * Dataset information struct
 * @param seed: dataset init seed
 * @param size: dataset size: 2^size
 * @param numClusters: amount of clusters
 * @param dataDim: amount of dimensions
*/
class Dataset {
    private:
        int seed; // generator seed
        int size; // dataset size 2^{size}
        unsigned int numClusters; // dataset cluster count
        unsigned int dataDim; // data dimension
        uint64_t numDataPoints; // amount of datapoints
        std::string outdir; // output directory
        std::default_random_engine generator; // generator used for generating points
        std::mt19937 mean_generator; // generator used for initial mean sets

        // Dataset vectors
        std::vector<std::vector<double>> clusterCenters;  // 2d vector of cluster centers
        std::vector<double> clusterStd;  // vector containing std dev of each cluster
        std::vector<uint64_t> clusterSize;  // vector containing amount of points in each cluster
        std::vector<uint32_t> membership;  // vector denoting intended cluster of datapoint at same index
        std::vector<std::vector<double>> datapoints;  // vector containing the datapoints
        std::vector<std::vector<uint64_t>> initialMeans;  // vector containing indexes of initial means
    
    public:
        // Constructor and Destructor
        Dataset(const int seed_, const int size_, const int numClusters_, const int dataDim_, const char *outdir_);
        ~Dataset();
        
        /**
         * Prints the dataset
        */
        void print() const;
        /**
         * Uniform randomly generate the means of the clusters in the interval [CENTRE_MIN, CENTRE_MAX]
        */
        void generate_cluster_means();
        /**
         * Uniform randomly generate the standard deviation of the clusters 
         * between 1/16 and 1/8 of the width of the interval [CENTRE_MIN, CENTRE_MAX]
        */
        void generate_cluster_std_devs();
        /**
         * Create a dataset using the provided dataset object and random generator 
        */
        void generate_dataset();
        /**
         * Upscale provided dataset from 2^x to 2^(x+y)
         * @param y: power increase of dataset size
        */
        void upscale_dataset(const unsigned int y);
        /**
         * Euclidean distance between two points.
         * Dimensions between points x and y need to be identical.
         * @param x: point 1
         * @param y: point 2
        */
        double distance(const std::vector<double> x, const std::vector<double> y) const;
        /**
         * Generate (set) of initial means to be used during execution.
         * This function uses the k-means++ method to generate initial means.
         * First centroid is randomly selected (uniformly).
         * In short, next centroid is the point with maximum distance from previously selected centroids.
         * This is different from the original k-means++ probability distribution, but achieves comparable results.
         * Selected centroids are used as initial means.
         * Initial means are saved in initial_means vector of dataset.
         * Based on implementation from https://www.geeksforgeeks.org/ml-k-means-algorithm/.
        */
        void generate_initial_means();
        /**
         * Write data points and their intended membership to corresponding files
        */
        int write_data() const;
        /**
         * Write generated cluster to file
        */
        int write_cluster_centers() const;
        /**
         * Write generated means to files
        */
        int write_initial_means() const;
        /**
         * Write dataset variables as metadata file
        */
        int write_metadata() const;
        /**
         * Calls all writing functions in one go
        */
        void write() const;
};

#endif /* __GENERATOR_H__ */
