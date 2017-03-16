#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "../common.h"

#include "../fastq_reader.h"
#include "../kmer_utils.h"
#include "../kmer_counter.h"
#include "../bloom_filter.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <map>

void deleteGeneratedFiles()
{
	// A lot of generated sequence file, lets remove it
	auto three_char_sequences = generateAllThreeCharSequenceCombinations();
	for (auto &file_name : three_char_sequences) 
	{
#ifdef USE_GZIP 
		remove((file_name + ".txt.gz").c_str());
#else
		remove((file_name + ".txt").c_str());
#endif
	}
}

TEST_CASE("test_fastq_reader", "[fastq_reader]") {
	
	FastqReader *reader = new FastqReader("1.fastq");
	auto kmer_files = reader->generateKmers(10);

	SECTION("number of returned kmer files") {
		REQUIRE(kmer_files.size() == 75);
	}

	SECTION("kmer file names contains specific names") {
		auto first_file_name = kmer_files.front()->filename;
		REQUIRE(first_file_name.find("AAAAA") != std::string::npos);

		auto last_file_name = kmer_files.back()->filename;
		REQUIRE(last_file_name.find("TTTTT") != std::string::npos);
	}
}

TEST_CASE("test_kmer_counter", "[kmer_counter]") {
	FastqReader *reader = new FastqReader("1.fastq");
	size_t kmersize = 30;
	size_t topcount = 25;
	auto kmer_files = reader->generateKmers(kmersize);
	auto first_kmer_file = kmer_files.front();
	auto kmer_counter = new TopKmerCounter(first_kmer_file, kmersize, topcount);
	auto map = kmer_counter->findTopKmers();

	SECTION("top kmer list size") {
		REQUIRE(map.size() == 4);
	}
	std::multimap<size_t, uint64_t>::iterator it = map.begin();
	SECTION("Top kmer count is 2") {
		REQUIRE(it->first == 2);
	}
	
	delete kmer_counter;
	delete reader;
}

TEST_CASE("delete_generated_files", "[delete]") {
	deleteGeneratedFiles();
}


/*
 *
void Test::test_fastaq_reader()
{
	FastqReader *reader = new FastqReader("1.fastq");
	std::string line;
	int count = 0;
	auto sequences = reader->generateSequences(10);
}

void printKmers(const std::vector<std::string> kmers)
{
	for (auto kmer : kmers)
	{
		std::cout << kmer << std::endl;
	}
}

void Test::test_find_kmers()
{
	std::string sequence = "TGATGGAACAGTGAAAGA";
	std::vector<std::string> kmers = findKmers(sequence, 10);
	printKmers(kmers);
	kmers.clear();
	kmers = findKmers(sequence, sequence.size());
	printKmers(kmers);
	kmers.clear();
	kmers = findKmers(sequence, 1);
	printKmers(kmers);
	kmers.clear();
	kmers = findKmers(sequence, sequence.size() - 1);
	printKmers(kmers);
	kmers.clear();

	sequence = "TGATGGAACAGTGAAAGATGAGACAAGCCCCGTGGAGGAGTGTTTTTTTAGTCAAAGTTCAAACTNATATCAGTGTCATACCATCACTNN";
	kmers = findKmers(sequence, 10);
	printKmers(kmers);
	kmers.clear();
}
	
void Test::test_find_top_kmers()
{
}
	
void Test::test_murmur_hash()
{
	
}

void Test::test_bloom_filter()
{
	std::string s1 = "ankara";
	std::string s2 = "istanbul";
	std::string s3 = "izmir";
	std::string s4 = "frankfurt";

	BloomFilter *bm = new BloomFilter(1000000, 10);
	bm->add(s1);
	bm->add(s2);
	bm->add(s3);
	bool b = bm->contains(s1);
	std::cout << s1 << ": " << b << std::endl;
	b = bm->contains(s2);
	std::cout << s2 << ": " << b << std::endl;
	b = bm->contains(s3);
	std::cout << s3 << ": " << b << std::endl;
	b = bm->contains(s4);
	std::cout << s4 << ": " << b << std::endl;
}
	
void Test::test_insert_kmer_to_map()
{
	
}
	
void Test::test_convert_to_kmer_long()
{
	std::string kmer = "AAAAAACCCAAACGGC";
	uint64_t expected = 1478194599297024;
	uint64_t kmer_long = convertToKmerLong(kmer);
	if (kmer_long != expected)
	{
		std::cout << "Error: " << kmer_long << " != " << expected << std::endl;
	}
	kmer = "A";
	expected = 0;
	kmer_long = convertToKmerLong(kmer);
	if (kmer_long != expected)
	{
		std::cout << "Error: " << kmer_long << " != " << expected << std::endl;
	}
	kmer = "AA";
	expected = 0;
	kmer_long = convertToKmerLong(kmer);
	if (kmer_long != expected)
	{
		std::cout << "Error: " << kmer_long << " != " << expected << std::endl;
	}
	kmer = "C";
	expected = 1;
	kmer_long = convertToKmerLong(kmer);
	if (kmer_long != expected)
	{
		std::cout << "Error: " << kmer_long << " != " << expected << std::endl;
	}
	kmer = "T";
	expected = 3;
	kmer_long = convertToKmerLong(kmer);
	if (kmer_long != expected)
	{
		std::cout << "Error: " << kmer_long << " != " << expected << std::endl;
	}
	kmer = "AT";
	expected = 3;
	kmer_long = convertToKmerLong(kmer);
	if (kmer_long != expected)
	{
		std::cout << "Error: " << kmer_long << " != " << expected << std::endl;
	}
}

void Test::test_convert_to_kmer_string()
{
	convertToKmerString(1478194599297024);
}
 */