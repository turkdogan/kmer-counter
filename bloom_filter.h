#ifndef BLOOM_FILTER_H
#define BLOOM_FILTER_H

#include <string>
#include <vector>

class BloomFilter {

public:
	BloomFilter(uint64_t size, uint8_t numHashes, uint8_t);
	~BloomFilter();
	/*
	 *
	void add(const std::string &) ;
	bool contains(const std::string &) const;
	 */
	size_t add(uint64_t);
	bool contains(uint64_t) const;
	size_t count(uint64_t);

private:
	/*
	 * There are some complaints about the vector<bool> on the internet.
	 * It is claimed that it does not guarantee contiguous memory.
	 http://stackoverflow.com/questions/19562667/how-to-zero-a-vectorbool
	 */
	std::vector<bool> m_bits;

	std::vector<uint8_t> m_counter;

	uint8_t m_max_count_value;

	uint8_t m_numHashes;

	uint64_t *m_buffer;
};
#endif
