CXX=g++
CXXFLAGS=-g -std=c++11 -Wall -pedantic
LIBS=../fastq_reader.o ../kmer_utils.o ../bloom_filter.o ../MurmurHash3.o ../kmer_counter.o -lz

SRC=tests.cpp
OBJ=$(SRC:%.cpp=%.o)

all: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $(LIBS) -o tests

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -f *.o	
	rm -f *.o	
