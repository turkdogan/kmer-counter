#include "bloom_filter.h"
#include "MurmurHash3.h"
#include <array>
#include <iostream>

inline uint64_t nthHash(uint8_t n,
	uint64_t hashA,
	uint64_t hashB,
	uint64_t filterSize) {
	return (hashA + n * hashB) % filterSize;
}

BloomFilter::BloomFilter(uint64_t size, uint8_t numHashes)
	: m_bits(size),
	m_numHashes(numHashes)
{
	// std::fill(m_counters.begin(), m_counters.end(), 0);
}


std::array<uint64_t, 2> hash(const uint8_t *data, std::size_t len) {
	std::array<uint64_t, 2> hash_value;
	MurmurHash3_x64_128(data, (int)len, 0, hash_value.data());

	return hash_value;
}

void BloomFilter::increment(size_t index)
{
	unsigned int i = m_counters[index];
	unsigned int mask = 1;
	while (i & mask)
	{
		i &= ~mask;
		mask <<= 1;
	}
	i |= mask;
	m_counters[index] = i;
}

void BloomFilter::add(const std::string &data) {
	auto hash_values = hash((uint8_t *)data.c_str(), data.size());

	for (int n = 0; n < m_numHashes; n++) {
		m_bits[nthHash(n, hash_values[0], hash_values[1], m_bits.size())] = true;
	//	increment(nthHash(n, hash_values[0], hash_values[1], m_counters.size()));
	}
}

bool BloomFilter::contains(const std::string &data) const {
	auto hash_values = hash((uint8_t *)data.c_str(), data.size());

	for (int n = 0; n < m_numHashes; n++) {
		if (!m_bits[nthHash(n, hash_values[0], hash_values[1], m_bits.size())]) {
			return false;
		}
	}
	return true;
}

BloomFilter::~BloomFilter()
{
	m_bits.clear();
}


