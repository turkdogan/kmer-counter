#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include "common.h"

#include <fstream>
#include <vector>
#include <map>

#ifdef USE_GZIP 
#include <zlib.h>
#endif

struct KmerFile
{
	std::string filename;
	size_t kmercount;
};

class FastqReader {

public:
	FastqReader(const char *filename);
	FastqReader(const char *filename, bool);
	std::vector<KmerFile *> generateKmers(size_t kmersize);

	bool readNextSequence(std::string &line);
private:

	std::ifstream m_file;


#ifdef USE_GZIP 
	std::map<std::string, gzFile> filename_map;
#else
	std::map<std::string, std::ofstream *> filename_map;
#endif

	std::map<std::string, KmerFile *> filename_fastq_map;
};
#endif
