/**
 * MPI kmeans dataset generator main file. 
 * This generator constructs a dataset of a set size,
 * including initial means to be supplied to the algorithm.
 * 
 * Author: Dennis Buurman
 * Based on work of Anne Hommelberg 'generateData.cc'
*/

#include "generator.h"

#include <iostream>
#include <random>
#include <sys/stat.h>

bool DEBUG = false;

/**
 * Checks if path is a directory
 * @param path: path to directory
*/
bool isdir(const char *path)
{
  struct stat statbuf;
  if (stat(path, &statbuf) != 0)
    return false;

  return S_ISDIR(statbuf.st_mode) != 0;
}

/**
 * Function for testing the dataset generations.
 * Running this function 2 times should result in identical datasets.
 * Only the initial mean sets should be different.
 * @param outdir: directory used for dataset initialization
*/
void test(const char *outdir) {
    int seed = 971, size = 8, clusters = 4, dim = 4, upscale = 1, mean_sets = 3;

    Dataset *d = new Dataset(seed, size, clusters, dim, outdir);
    d->generate_cluster_means();
    d->generate_cluster_std_devs();
    d->generate_dataset();
    d->upscale_dataset(upscale);
    for (int i = 0; i < mean_sets; i++) {
        d->generate_initial_means();
    }

    d->print();
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

    // DEBUGGING ONLY
    if (DEBUG) {
        std::cout << "WARNING: DEBUG MODE ACTIVE!" << std::endl;
        test(outdir);
        return 0;
    }

    std::default_random_engine generator(seed);
    Dataset *dataset = new Dataset(seed, MINIMUM_SIZE, numClusters, dataDim, outdir);

    std::cout << "Generating size " << size << " dataset with data of dimension " << dataDim << " in " << numClusters << " clusters using random seed " << seed << ". Writing output to " << outdir << std::endl;
    
    // Generate cluster means and standard deviations
    dataset->generate_cluster_means();
    dataset->generate_cluster_std_devs();
    // Generate default dataset
    dataset->generate_dataset();
    // Generate sets of initial means that can be supplied to algorithms
    std::cout << "Generating " << MEAN_SETS << " sets of initial means..." << std::endl;
    for (int i = 0; i < MEAN_SETS; i++) {
        dataset->generate_initial_means();
    }
    // Resize dataset by duplication
    if (size_difference > 0) {
        dataset->upscale_dataset(size_difference);
    } else if (size_difference < 0) {
        std::cerr << "ERROR: requested size " << size << " smaller than minimum size " << MINIMUM_SIZE << std::endl;
    }
    // Write dataset to files
    std::cout << "Writing output to files..." << std::endl;
    dataset->write();

    delete dataset;

    return 0;
}