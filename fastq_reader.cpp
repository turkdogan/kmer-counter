#include "fastq_reader.h"
#include "kmer_utils.h"

#include <string>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

FastqReader::FastqReader(const char* filename)
	: m_file(filename)
{
	if (!m_file.is_open())
	{
		throw std::invalid_argument("File could not be opened: " + std::string(filename));
	}
}

size_t FastqReader::getApproximateKmercount(size_t kmersize)
{
	// we are interesting with each 1 of 4 lines
	float kmer_line_percentage = 0.25;
	std::string line;
	readNextSequence(line);
	size_t line_count = std::count(std::istreambuf_iterator<char>(m_file),
		std::istreambuf_iterator<char>(), '\n');
	size_t kmer_line_count = line_count * kmer_line_percentage;
	size_t approximate_kmer_per_line = (line.size() - kmersize);
	return approximate_kmer_per_line * kmer_line_count;
}

bool FastqReader::readNextSequence(std::string &line)
{
	/*
	 * Assumed that, the structure of file is OK
	 * Therefore no defensive code written
	 */
	m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	auto b = std::getline(m_file, line) ? true : false;
	m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	return b;
}

void FastqReader::reload()
{
	m_file.clear();
	m_file.seekg(0, m_file.beg);
}

