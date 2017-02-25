#ifndef KMER_UTILS_H
#define KMER_UTILS_H

#include "common.h"

#include <vector>
#include <string>
#include <map>
	
#ifdef USE_GZIP 
#include <zlib.h>
#include <errno.h>
#include <string.h>
#endif

std::vector<std::uint64_t> findKmers(const std::string &, size_t kmersize);

void insert_kmer_to_map(std::string &kmer, std::map<uint64_t, int> &map);

uint64_t encodeSequence(std::string &kmer);

char bitstoc(unsigned int twobit);

/*
 * Kmer is encoded as unsigned long long
 */
void decodeSequence(uint64_t seq, char *buf, size_t size);
std::string decodeSequence(uint64_t seq, size_t size);

std::vector<std::string> generateAllThreeCharSequenceCombinations();
	
	#ifdef USE_GZIP 

		gzFile openGzipFile(const char *filename, const char *readWrite);
		void closeGzipFile(gzFile file);
		size_t readcompressedFile(gzFile file, void *out, size_t len);
		void writeCompressedFile(gzFile file, const void *buffer, size_t len);

	#endif


#endif