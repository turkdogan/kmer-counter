#include "kmer_counter.h"
#include "kmer_utils.h"
#include <iostream>

#define SMALL_FILE_KMER_COUNT 10000

namespace
{
	template <typename A, typename B>
	std::multimap<B, A> flip_map(std::map<A, B> & src)
	{
		std::multimap<B, A> destination;

		for (auto it = src.begin(); it != src.end(); ++it)
			destination.insert(std::pair<B, A>(it->second, it->first));

		return destination;
	}
}

TopKmerCounter::TopKmerCounter(FastqReader *reader, size_t kmersize = 30, size_t topcount = 25, size_t kmignore = 0)
	: m_file(reader), m_kmersize(kmersize), m_topcount(topcount)
{
	m_approx_kmer_count = reader->getApproximateKmercount(kmersize);
	reader->reload();
	if (kmignore == 0)
	{
		size_t small_file_kmer_size = SMALL_FILE_KMER_COUNT;
		if (m_approx_kmer_count < small_file_kmer_size)
		{
			m_ignored_kmer_count = 1;
		}
		else if (m_approx_kmer_count < small_file_kmer_size * 10)
		{
			m_ignored_kmer_count = 4;
		}
		else if (m_approx_kmer_count < small_file_kmer_size * 100 * 100)
		{
			m_ignored_kmer_count = 10;
		}
		else if (m_approx_kmer_count < small_file_kmer_size * 1000 * 100)
		{
			m_ignored_kmer_count = 50;
		}
		else
		{
			m_ignored_kmer_count = 255;
		}
	}
	else
	{
		m_ignored_kmer_count = kmignore;
	}
}

/*
 * Returns topcount k-mers in the given sequence file
 * The unique kmers are ignored.
 */
std::multimap<size_t, uint64_t> TopKmerCounter::findTopKmers()
{
	processMapWithBloomFilter();
	postprocess();
	return m_final_map;
}

void TopKmerCounter::processMapWithBloomFilter()
{
	std::cout << "ignore: " << m_ignored_kmer_count << std::endl;
	m_bloom_filter = new BloomFilter(m_approx_kmer_count * 5, 5, m_ignored_kmer_count);
	std::string line;

	while (m_file->readNextSequence(line))
	{
		auto buffer = findKmers(line, m_kmersize);
		for (size_t i = 0; i < buffer.size(); i++)
		{
			uint64_t kmer = buffer[i];
			size_t response = m_bloom_filter->add(kmer);
			if (response != 0)
			{
				if (m_kmer_to_count_map.count(kmer) == 0)
				{
					m_kmer_to_count_map[kmer] = response + 1;
				}
				else
				{
					m_kmer_to_count_map[kmer]++;
				}
			}
		}
		buffer.clear();
	}
}

/*
 * Sort map
 */
void TopKmerCounter::postprocess()
{
	// sort map by count
	auto sorted = flip_map(m_kmer_to_count_map);

	auto top_counter = m_topcount;

	// map is sorted by ascending order
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
		m_final_map.insert(std::make_pair(iter->first, iter->second));
		++iter;

		// # of kmemrs can be more than m_topcount
		// for example the following map is legal for m_topcount = 2
		// "foo1" = 10, "foo2" = 10, "foo3" = 10, "foo4" = 3
		if (iter != sorted.rend() && last_max_count > iter->first)
		{
			last_max_count = iter->first;
			top_counter--;
		}
	}
	m_kmer_to_count_map.clear();
}

TopKmerCounter::~TopKmerCounter()
{
	delete m_bloom_filter;
}
