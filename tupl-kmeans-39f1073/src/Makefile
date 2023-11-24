CXX = mpic++
CXXFLAGS = -Wall -Wextra -O3 -lrt -std=c++0x -g

TARGETS = \
	MPI_Kmeans		\
	MPI_Kmeans_incremental	\
	MPI_Kmeans_localized	\
	MPI_Kmeans_incremental_localized

HEADERS = mpi-utils.h common.h algorithm.h

OBJS = common.o


all:	$(TARGETS)

MPI_Kmeans:	MPI_Kmeans.cpp $(OBJS) $(HEADERS)
		$(CXX) $(CXXFLAGS) -o $@ $^

MPI_Kmeans_localized:	MPI_Kmeans_localized.cpp $(OBJS) $(HEADERS)
			$(CXX) $(CXXFLAGS) -o $@ $^


MPI_Kmeans_incremental:	MPI_Kmeans_incremental.cpp $(OBJS) $(HEADERS)
			$(CXX) $(CXXFLAGS) -o $@ $^

MPI_Kmeans_incremental_localized:	MPI_Kmeans_incremental_localized.cpp $(OBJS) $(HEADERS)
					$(CXX) $(CXXFLAGS) -o $@ $^


%.o:	%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(TARGETS)
	rm -f $(OBJS)
