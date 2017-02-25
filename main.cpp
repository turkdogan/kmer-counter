#include "common.h"
#include "kmer_counter.h"
#include "kmer_utils.h"
#include <iostream>
#include <ctime>
#include <algorithm>
#include <stdio.h>
#include <mutex>
#include <future>
#include <sstream>


using namespace std;

std::multimap<size_t, uint64_t> global_result;

void getTopKmers(KmerFile *kmer_file, size_t kmersize, size_t topcount)
{
	auto kmer_counter = new TopKmerCounter(kmer_file, kmersize, topcount);
	auto map = kmer_counter->findTopKmers();
	global_result.insert(map.begin(), map.end());
	delete kmer_counter;
	map.clear();
#ifdef USE_GZIP 
		// std::remove((kmer_file->filename + ".gz").c_str());
#else
		std::remove((kmer_file->filename).c_str());
#endif
}

void runMultiThread(const char *filename, size_t kmersize, size_t topcount, const char *outfilename)
{
	
}

void runSingleThread(const char *filename, size_t kmersize, size_t topcount, const char *outfilename)
{
	
	auto fastq_reader = new FastqReader(filename);
	auto kmer_files = fastq_reader->generateKmers(kmersize);
	delete fastq_reader;

	for (auto &kmer_file : kmer_files)
	{
		getTopKmers(kmer_file, kmersize, topcount);
	}

	ofstream out(outfilename);

	auto iter = global_result.rbegin();

	// find m_topcount elements (including duplicates)
	while (iter != global_result.rend() && topcount > 0)
	{
		auto last_max_count = iter->first;
		std::cout << iter->first << " : " << decodeSequence(iter->second, kmersize)<< std::endl;
		out << iter->first << " : " << decodeSequence(iter->second, kmersize) << std::endl;

		++iter;

		if (iter->first < last_max_count)
		{
			last_max_count = iter->first;
			topcount--;
		}
	}
	std::cout << "Results are in the " << outfilename << std::endl;

	out.close();

	global_result.clear();
}

void runApplication(const char *filename, size_t kmersize, size_t topcount, const char *outfilename)
{
	std::cout << "filename: " << filename << endl;
	std::cout << "kmer size: " << kmersize << endl;
	std::cout << "top count: " << topcount << endl;

	runSingleThread(filename, kmersize, topcount, outfilename);
	// runMultiThread(filename, kmersize, topcount, outfilename);

}

// http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
	auto itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
		return *itr;
	}
	return nullptr;
}

void printOptions()
{
	std::cout << "parameters should be: --filename ${filename} --kmersize ${size} --topcount {count}" << std::endl;
}

/*
 * 		gzFile openGzipFile(const char *filename, const char *readWrite);
		void closeGzipFile(gzFile file);
		size_t readcompressedFile(gzFile file, void *out, size_t len);
		void writeCompressedFile(gzFile file, const void *buffer, size_t len);

 */

void foo()
{
	uint64_t kmer = 0;
	size_t len = 0;

	std::string s = "GTAAACCCGTGTGTACGTAAACCCGTGTGT";
	uint64_t encoded = encodeSequence(s);
	uint64_t encoded_read = encodeSequence(s);

#ifdef USE_GZIP 
	gzFile file = openGzipFile("foo.gz", "w");
	writeCompressedFile(file, &encoded, sizeof(encoded));
	closeGzipFile(file);

	file = openGzipFile("foo.gz", "r");
	while ((len = readcompressedFile(file, &encoded_read, sizeof(encoded_read))) > 0)
	{
		std::cout << encoded << " " << encoded_read << " " << std::endl;
		std::cout <<  decodeSequence(encoded_read, s.length()) << std::endl;
		std::cout <<  s << std::endl;
	}
	closeGzipFile(file);
#endif

}

int main(int argc, char **argv) 
{
	auto filename = getCmdOption(argv, argv + argc, "--filename");
	auto kmersize = getCmdOption(argv, argv + argc, "--kmersize");
	auto topcount = getCmdOption(argv, argv + argc, "--topcount");

	if (!filename || !kmersize || !topcount)
	{
		printOptions();
		return 0;
	}
	char *outfilename = getCmdOption(argv, argv + argc, "--outfilename");
	if (!outfilename)
	{
		outfilename = (char *)"out.txt";
	}

	try {
		auto begin_time = clock();
		runApplication(filename, std::stoi(kmersize), std::stoi(topcount), outfilename);
		std::cout << "Total duration: " << float(clock() - begin_time) / (CLOCKS_PER_SEC * 60) << " minutes" << std::endl;
	}
	catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		// TODO clear objects here (files etc.)
		throw;
	}
}