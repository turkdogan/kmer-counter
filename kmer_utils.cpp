#include "kmer_utils.h"
#include <iostream>
#include <sstream>

#ifdef USE_GZIP 

gzFile openGzipFile(const char *filename, const char *readWrite)
{
	gzFile file = gzopen(filename, readWrite);
	if (!file) {
		fprintf(stderr, "gzopen of '%s' failed: %s.\n", filename,
			strerror(errno));
		exit(EXIT_FAILURE);
	}
	return file;
}

void closeGzipFile(gzFile file)
{
	gzclose(file);
}

size_t readcompressedFile(gzFile file, void *out, size_t len)
{
	int err;
	size_t bytes_read;
	bytes_read = gzread(file, out, len);
	// printf ("%s", buffer);
	if (bytes_read < len) {
		if (!gzeof(file)) {
			const char * error_string;
			error_string = gzerror(file, &err);
			if (err) {
				fprintf(stderr, "Error: %s.\n", error_string);
				exit(EXIT_FAILURE);
			}
		}
	}
	return bytes_read;
}

void writeCompressedFile(gzFile file, const void *buffer, size_t len)
{
	int bytes_written = gzwrite(file, buffer, len);
	if (bytes_written == 0)
	{
		int err_no = 0;
		fprintf(stderr, "Error during compression: %s", gzerror(file, &err_no));
		gzclose(file);
	}
}

#endif

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
		//kmers.push_back(sequence.substr(i, kmersize));
	}
	return kmers;
}

void insert_kmer_to_map(std::string &kmer, std::map<uint64_t, int> &map)
{
	// TODO check it
	// TODO write test
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

std::vector<std::string> generateAllThreeCharSequenceCombinations()
{
	std::vector<std::string> sequence_combinations;
	std::string genomChars("ACGT");
	for (auto c1 : genomChars)
	{
		std::string s1(1, c1);
		for (auto c2 : genomChars)
		{
			std::string s2(1, c2);
			for (auto c3 : genomChars)
			{
				std::string s3(1, c3);
				for (auto c4 : genomChars)
				{
					std::string s4(1, c4);
					for (auto c5: genomChars)
					{
						std::string s5(1, c5);
						sequence_combinations.push_back(s1 + s2 + s3 + s4 + s5);
					}
				}
			}
		}
	}
	return sequence_combinations;
}
