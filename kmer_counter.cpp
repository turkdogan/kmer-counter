#include "kmer_counter.h"
#include "kmer_utils.h"
#include <iostream>

// use bloom filter if more than BLOOM_FILTER_THRESHOLD
#define BLOOM_FILTER_THRESHOLD 1000

namespace
{
	template <typename A, typename B>
	std::multimap<B, A> flip_map(std::map<A, B> & src) {

		std::multimap<B, A> destination;

		for (auto it = src.begin(); it != src.end(); ++it)
			destination.insert(std::pair<B, A>(it->second, it->first));

		return destination;
	}
}

TopKmerCounter::TopKmerCounter(KmerFile *kmerfile, size_t kmersize = 30, size_t topcount = 25)
#ifdef USE_GZIP 
	: m_kmercount_in_sequence_file(kmerfile->kmercount), m_kmersize(kmersize), m_topcount(topcount)
#else
	: m_file(kmerfile->filename), m_kmercount_in_sequence_file(kmerfile->kmercount), m_kmersize(kmersize), m_topcount(topcount)
#endif
{
#ifdef USE_GZIP 
	std::string compressed_file(kmerfile->filename + ".gz");
	m_file = openGzipFile(compressed_file.c_str(), "r");
#endif
}

/*
 * Returns topcount k-mers in the given sequence file
 * The unique kmers are ignored. 
 */
std::multimap<size_t, std::string> TopKmerCounter::findTopKmers()
{
	if (m_kmercount_in_sequence_file < BLOOM_FILTER_THRESHOLD)
	{
		processMapWithoutBloomFilter();
	}
	else
	{
		processMapWithBloomFilter();
	}
	mergeBucketKmers();
	return m_final_map;
}

void TopKmerCounter::processMapWithBloomFilter()
{
	m_bloom_filter = new BloomFilter(m_kmercount_in_sequence_file * 100, 4);
	std::string kmer;

#ifdef USE_GZIP 
	size_t len;
	char *buffer = new char[m_kmersize + 1];
	while ((len = readcompressedFile(m_file, buffer, m_kmersize, true)) > 0)
	{
		kmer = std::string(buffer);
#else
	while (std::getline(m_file, kmer))
	{
#endif
		// A couple of reason to "trust" bloom filter:
		// 1 - False positive is ratio is very low, 
		// because in every BUCKET_SIZE, bloom filter is cleared
		// 2 - False positive is ignorable because,
		// we are finding the top "n" count of kmers.
		// the false positive k-mers cannot be in the top k-mers
		if (m_bloom_filter->contains(kmer))
		{
			// if the bloom filter already contains
			// add it to the counter hash table
			// m_kmer_to_count_map[kmer]++;
		} 
		else
		{
			// the second encountered k-mer is inserted into the bloom filter
			// it is not inserted to the hash map. As most of the kmers in a 
			// genome is unique, which speeds up the code. When calculating 
			// the final counts, we add +1 to the result for this kmer.
			m_bloom_filter->add(kmer);
		}
	}
#ifdef USE_GZIP 
	closeGzipFile(m_file);
#else
	m_file.close();
#endif
	delete m_bloom_filter;
}

void TopKmerCounter::processMapWithoutBloomFilter()
{
	std::string kmer;
#ifdef USE_GZIP 
	size_t len;
	char *buffer = new char[m_kmersize + 1];
	while ((len = readcompressedFile(m_file, buffer, m_kmersize, true)) > 0)
	{
		kmer = std::string(buffer);
#else
	while (std::getline(m_file, kmer))
	{
#endif
		m_kmer_to_count_map[kmer]++;
	}
}

/*
 * On every BUKET_SIZE k-mer top m_topcount kmers are kept
 * on a temp map. This map may contain more than m_topcount kmers.
 * Because some k-mers count can be equal to each other. Both of them
 * must be kept
 */
void TopKmerCounter::mergeBucketKmers()
{
	// sort map by count
	auto sorted = flip_map(m_kmer_to_count_map);

	auto top_counter = m_topcount;

	// map is sorted by ascending order,
	// to iterate from end for descending count order
	auto iter = sorted.rbegin();

	// find m_topcount elements (including duplicates)
	while (iter != sorted.rend() && top_counter > 0)
	{
		auto last_max_count = iter->first;

		// may have duplicate kmer from the previous
		// bucket processings, need to merge those duplicates
		// +1 for the kmer in the bloom filter
		//char *decoded;
		//decodeSequence(iter->second, decoded, m_kmersize);
		m_final_map.insert(std::make_pair(iter->first + 1, iter->second));
		++iter;

		// # of kmemrs can be more than m_topcount
		// for example the following map is legal for m_topcount = 2
		// "foo1" = 10, "foo2" = 10, "foo3" = 10, "foo4" = 3
		if (last_max_count > iter->first)
		{
			last_max_count = iter->first;
			top_counter--;
		}
	}
	m_kmer_to_count_map.clear();
}

TopKmerCounter::~TopKmerCounter()
{
}