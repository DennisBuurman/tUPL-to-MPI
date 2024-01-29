CXX = g++
CXXFLAGS = -Wall -Wextra -g -O3 -std=c++0x


all:	generator

generator:	generator.cc
		$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f generator
