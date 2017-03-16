#include "kmer_counter.h"
#include "kmer_utils.h"
#include <iostream>
#include <ctime>
#include <sstream>
#include <algorithm>

using namespace std;

void runSingleThread(const char *filename, size_t kmersize, size_t topcount, const char *outfilename, size_t kmignore)
{
	auto fastq_reader = new FastqReader(filename);
	auto kmer_counter = new TopKmerCounter(fastq_reader, kmersize, topcount, kmignore);
	auto map = kmer_counter->findTopKmers();

	ofstream out(outfilename);

	auto iter = map.rbegin();

	// find m_topcount elements (including duplicates)
	while (iter != map.rend())
	{
		std::cout << iter->first << " : " << decodeSequence(iter->second, kmersize) << std::endl;
		out << iter->first << " : " << decodeSequence(iter->second, kmersize) << std::endl;
		++iter;
	}
	std::cout << "Results are in the " << outfilename << std::endl;

	out.close();

	delete kmer_counter;
	delete fastq_reader;
	map.clear();
}

void runApplication(const char *filename, size_t kmersize, size_t topcount, const char *outfilename, size_t kmignore)
{
	std::cout << "filename: " << filename << endl;
	std::cout << "kmer size: " << kmersize << endl;
	std::cout << "top count: " << topcount << endl;

	runSingleThread(filename, kmersize, topcount, outfilename, kmignore);
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
	std::cout << "parameters should be: --filename ${filename} --kmersize ${size} --topcount {count} --kmignore {size}" << std::endl;
}

int main(int argc, char **argv)
{
	auto filename = getCmdOption(argv, argv + argc, "--filename");
	auto kmersize = getCmdOption(argv, argv + argc, "--kmersize");
	auto topcount = getCmdOption(argv, argv + argc, "--topcount");
	auto kmignore = getCmdOption(argv, argv + argc, "--kmignore");

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
	size_t km_ignore = 0;
	if (kmignore)
	{
		km_ignore = std::stoi(kmignore);
	}

	try {
		auto begin_time = clock();
		runApplication(filename, std::stoi(kmersize), std::stoi(topcount), outfilename, km_ignore);
		std::cout << "Total duration: " << float(clock() - begin_time) / (CLOCKS_PER_SEC * 60) << " minutes" << std::endl;
	}
	catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		// TODO clear objects here (files etc.)
		throw;
	}
}
