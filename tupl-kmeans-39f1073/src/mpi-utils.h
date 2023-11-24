/*
 * MPI K-means implementation
 *
 * Author: Anne Hommelberg
 *
 * Code refactoring done by Kristian Rietveld, Leiden University.
 */

#ifndef __MPI_UTILS__
#define __MPI_UTILS__

#include <mpi.h>

#include <tuple>
#include <string>

#ifndef MPI_UINT64_T
#	define MPI_UINT64_T MPI_UNSIGNED_LONG
#endif

extern int mpi_rank, mpi_size;
extern char mpi_hostname[MPI_MAX_PROCESSOR_NAME];

extern MPI_Datatype MPI_DATAPOINT;


static inline std::string printPreamble(void)
{
  return std::string("[") + std::to_string(mpi_rank + 1) +
      std::string("/") + std::to_string(mpi_size) + std::string("]");
}

static inline std::pair<uint64_t, uint64_t> mpi_range(uint64_t size, int rank, int count)
{
  uint64_t part = size/count;
  uint64_t from = rank*part;
  uint64_t to = from+part;
  if (count == rank+1)
    to = size;

  return {from, to};
}

static inline std::pair<uint64_t, uint64_t> mpi_range(uint64_t size, int rank)
{
  return mpi_range(size, rank, mpi_size);
}

static inline std::pair<uint64_t, uint64_t> mpi_range(uint64_t size)
{
  return mpi_range(size, mpi_rank, mpi_size);
}

#endif /* __MPI_UTILS_H__ */
