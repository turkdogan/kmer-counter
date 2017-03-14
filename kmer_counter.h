#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include <map>
#include "fastq_reader.h"
#include "bloom_filter.h"
#include "common.h"
#include <vector>

class TopKmerCounter {

public:
	TopKmerCounter(KmerFile *kmerfile, size_t kmersize, size_t topcount);
	~TopKmerCounter();
	
	std::multimap<size_t, uint64_t> findTopKmers();


private:

#ifdef USE_GZIP 
	gzFile m_file;
#else
	std::ifstream m_file;
#endif
	size_t m_kmercount_in_sequence_file;

	size_t m_kmersize;
	size_t m_topcount;

	/*
	 * We keep an internal bloom filter 
	 * to check if the kmer already in the filter.
	 * If we have at least twice kmer in the filter
	 * then we store this value to the map
	 */
	BloomFilter *m_bloom_filter;

	std::map<uint64_t, size_t> m_kmer_to_count_map;
	std::multimap<size_t, uint64_t> m_final_map;

	void processMapWithBloomFilter();

	void mergeBucketKmers();
};

#endif
