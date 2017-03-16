#include "bloom_filter.h"
#include "MurmurHash3.h"
#include <array>
#include <iostream>

namespace
{
	void fillHashes(
		uint64_t *buffer,
		uint8_t size,
		uint64_t hashA,
		uint64_t hashB,
		uint64_t filterSize) {
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = (hashA + i * hashB) % filterSize;
		}
	}
}

BloomFilter::BloomFilter(uint64_t size, uint8_t numHashes, uint8_t max_count_value)
	: m_counter(size),
	m_max_count_value(max_count_value),
	m_numHashes(numHashes)
{
	m_buffer = new uint64_t[numHashes];
}

std::array<uint64_t, 2> hash(const uint8_t *data, std::size_t len) {
	std::array<uint64_t, 2> hash_value;
	MurmurHash3_x64_128(data, (int)len, 0, hash_value.data());

	return hash_value;
}

size_t BloomFilter::add(uint64_t data)
{
	auto hash_values = hash((uint8_t *)&data, sizeof(data));
	fillHashes(m_buffer, m_numHashes, hash_values[0], hash_values[1], m_counter.size());
	uint64_t min = m_max_count_value;
	uint64_t max = 0;
	for (size_t n = 0; n < m_numHashes; n++)
	{
		if (m_counter[m_buffer[n]] > max)
		{
			max = m_counter[m_buffer[n]];
		}
		if (m_counter[m_buffer[n]] < min)
		{
			min = m_counter[m_buffer[n]];
		}
	}
	if (max < m_max_count_value)
	{
		for (size_t n = 0; n < m_numHashes; n++)
		{
			size_t index = m_buffer[n];
			//m_bits[index] = true;
			m_counter[index]++;
		}
		return 0;
	}
	return min;
	/*
	 *
	std::cout << std::endl;
	for (int n = 0; n < m_numHashes; n++) {
		std::cout << m_bits[nthHash(n, hash_values[0], hash_values[1], m_bits.size())] << " ";
	}
	std::cout << std::endl;
	 */
}

bool BloomFilter::contains(uint64_t data) const
{
	auto hash_values = hash((uint8_t *)&data, sizeof(data));
	fillHashes(m_buffer, m_numHashes, hash_values[0], hash_values[1], m_counter.size());

	for (size_t n = 0; n < m_numHashes; n++)
	{
		if (m_counter[m_buffer[n]] == 0) {
			return false;
		}
	}
	return true;
}


/*
 *
void BloomFilter::add(const std::string &data) {
	auto hash_values = hash((uint8_t *)data.c_str(), data.size());

	for (int n = 0; n < m_numHashes; n++) {
		m_bits[nthHash(n, hash_values[0], hash_values[1], m_bits.size())] = true;
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
 */

BloomFilter::~BloomFilter()
{
	//m_bits.clear();
	m_counter.clear();
}


