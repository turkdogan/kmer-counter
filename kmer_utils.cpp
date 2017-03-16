#include "kmer_utils.h"
#include <iostream>
#include <sstream>

static inline uint8_t ctobits(char c)
{
	switch (c) {
	case 'A':
	case 'N': return 0;
	case 'C': return 1;
	case 'T': return 2;
	case 'G': return 3;
	default:
		fprintf(stderr, "Unexpected char: %c\n", c);
		exit(EXIT_FAILURE);
		break;
	}
}


char bitstoc(unsigned int twobit)
{
	static char seq_tab[] = { 'A', 'C', 'T', 'G' };
	return seq_tab[twobit];
}

/*
* s = std::string("ATCG")
* dnsseq_encode(s.c_str(), s.length())
*/
uint64_t dnaseq_encode(const char *s, size_t len)
{
	const char *cp = s, *end = s + len;
	unsigned int nr_shift = 0;
	uint64_t ret = 0;

	for (; cp < end; cp++, nr_shift += 2) {
		uint64_t v = ctobits(*cp);
		ret |= v << nr_shift;
	}
	return ret;
}

std::string decodeSequence(uint64_t seq, size_t size)
{
	std::stringstream out;
	size_t i = 0;
	while (i < size)
	{
		out << bitstoc(seq & 3U);
		seq = seq >> 2;
		i++;
	}
	return out.str();
}

void decodeSequence(uint64_t seq, char *buf, size_t size)
{
	const char *end = buf + size;

	for (; buf < end; buf++) {
		*buf = bitstoc(seq & 3U);
		seq = seq >> 2;
	}
	*buf = '\0';
}

size_t findKmers(const std::string &sequence, size_t kmersize, std::vector<uint64_t> &buffer)
{
	if (kmersize < 1 || kmersize > sequence.size())
	{
		throw std::invalid_argument(
			"Recieved invalid sequence or kmer_size: " + sequence +
			" kmer_size: " + std::to_string(kmersize));
	}
	std::vector<std::uint64_t> kmers;
	size_t size = sequence.size();
	size_t size_response = size - kmersize + 1;
	if (buffer.size() <= size_response)
	{
		buffer.resize(size_response);
	}
	for (size_t i = 0; i <= size - kmersize; i++)
	{
		std::string kmer = sequence.substr(i, kmersize);
		buffer[i] = encodeSequence(kmer);
	}
	return size_response;
}

std::vector<std::uint64_t> findKmers(const std::string &sequence, size_t kmersize)
{
	if (kmersize < 1 || kmersize > sequence.size())
	{
		throw std::invalid_argument(
			"Recieved invalid sequence or kmer_size: " + sequence +
			" kmer_size: " + std::to_string(kmersize));
	}
	std::vector<std::uint64_t> kmers;
	size_t size = sequence.size();
	for (size_t i = 0; i <= size - kmersize; i++)
	{
		std::string kmer = sequence.substr(i, kmersize);
		uint64_t encoded = encodeSequence(kmer);
		kmers.push_back(encoded);
	}
	return kmers;
}

void insert_kmer_to_map(std::string &kmer, std::map<uint64_t, int> &map)
{
	uint64_t key = encodeSequence(kmer);
	map[key]++;
}

uint64_t encodeSequence(std::string &kmer)
{
	const char *cp = kmer.c_str();
	const char *end = cp + kmer.size();
	unsigned int nr_shift = 0;
	uint64_t ret = 0;

	for (; cp < end; cp++, nr_shift += 2) {
		uint64_t v = ctobits(*cp);
		ret |= v << nr_shift;
	}
	return ret;
}
