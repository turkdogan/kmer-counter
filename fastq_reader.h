#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <fstream>
#include <vector>
#include <map>

class FastqReader {

public:
	FastqReader(const char *filename);

	bool readNextSequence(std::string &line);
	size_t getApproximateKmercount(size_t kmersize);

	void reload();
private:

	std::ifstream m_file;
};
#endif
