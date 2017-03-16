#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include <map>
#include "fastq_reader.h"
#include "bloom_filter.h"
#include <vector>

class TopKmerCounter {

public:
	TopKmerCounter(FastqReader *reader, size_t kmersize, size_t topcount, size_t kmignore);
	~TopKmerCounter();

	std::multimap<size_t, uint64_t> findTopKmers();


private:

	FastqReader *m_file;

	size_t m_kmersize;
	size_t m_topcount;

	size_t m_approx_kmer_count;

	size_t m_ignored_kmer_count;

	/*
	 * We keep an internal bloom filter
	 * to check if the kmer already in the filter.
	 * If we have at least twice kmer in the filter
	 * then we store this value to the map
	 */
	BloomFilter *m_bloom_filter;

	std::map<uint64_t, size_t> m_kmer_to_count_map;
	std::multimap<size_t, uint64_t> m_final_map;

	void postprocess();
	void mergeFilterWithMap();
	void processMapWithBloomFilter();
};

#endif
