CXX = g++
CXXFLAGS = -c -Wall -Wextra -g -O3 -std=c++0x
LDFLAGS =
SOURCES = main.cc generator.cc
OBJECTS = $(SOURCES:.cc=.o)
EXECUTABLE = generator

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

.cc.o:
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -rf *o $(EXECUTABLE)
