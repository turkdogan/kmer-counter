#ifndef BLOOM_FILTER_H
#define BLOOM_FILTER_H

#include <string>
#include <vector>
#include <bitset>

class BloomFilter {

public:
	BloomFilter(uint64_t size, uint8_t numHashes);
	~BloomFilter();
	void add(const std::string &) ;
	bool contains(const std::string &) const;
	void add(uint64_t) ;
	bool contains(uint64_t) const;

private:
	std::vector<unsigned short> m_counters;
	/*
	 * There is some complaints about the vector<bool> on the internet.
	 * It is claimed that it does not guarantee contiguous memory.
	 http://stackoverflow.com/questions/19562667/how-to-zero-a-vectorbool
	 */
	std::vector<bool> m_bits;
	
	uint8_t m_numHashes;

	void increment(size_t index);

};
#endif
