#ifndef KMER_UTILS_H
#define KMER_UTILS_H

#include <vector>
#include <string>
#include <map>
	
std::vector<std::uint64_t> findKmers(const std::string &, size_t kmersize);

size_t findKmers(const std::string &sequence, size_t kmersize, std::vector<uint64_t> &buffer);

void insert_kmer_to_map(std::string &kmer, std::map<uint64_t, int> &map);

uint64_t encodeSequence(std::string &kmer);

char bitstoc(unsigned int twobit);

/*
 * Kmer is encoded as unsigned long long
 */
void decodeSequence(uint64_t seq, char *buf, size_t size);
std::string decodeSequence(uint64_t seq, size_t size);

#endif
