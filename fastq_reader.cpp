#include "fastq_reader.h"
#include "kmer_utils.h"

#include <string>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <sstream>

FastqReader::FastqReader(const char* filename)
	: m_file(filename)
{
	if (!m_file.is_open())
	{
		throw std::invalid_argument("File could not be opened: " + std::string(filename));
	}
	auto three_char_sequences = generateAllThreeCharSequenceCombinations();
#ifdef USE_GZIP 
	for (auto &seq : three_char_sequences)
	{
		filename_map[seq] = openGzipFile(std::string("temp/" + seq + ".txt.gz").c_str(), "wb");
	}
#else
	for (auto &seq : three_char_sequences)
	{
		filename_map[seq] = new std::ofstream("temp/" + seq + ".txt");
	}
#endif
}

void FastqReader::insert_sequence(std::map<std::string, std::ostringstream> &buffer_map)
{
	for (auto &pair : buffer_map)
	{
		if (!pair.first.empty())
		{
#ifdef USE_GZIP 
			std::string s = pair.second.str();
			writeCompressedFile(filename_map[pair.first], s.c_str(), s.size());
#else
			*filename_map[pair.first] << pair.second.str();
#endif
			if (filename_fastq_map[pair.first])
			{
				filename_fastq_map[pair.first]->kmercount++;
			} 
			else
			{
				auto fastqFile = new KmerFile();
				fastqFile->filename = std::string("temp/" + pair.first + ".txt").c_str();
				fastqFile->kmercount = 1;
				filename_fastq_map[pair.first] = fastqFile;
			}
		}
	}
}

static std::string getfirstFiveChars(uint64_t kmer) 
{
	unsigned int c1 = kmer & 0x3;
	unsigned int c2 = (kmer >> 2) & 0x3;
	unsigned int c3 = (kmer >> 4) & 0x3;
	unsigned int c4 = (kmer >> 6) & 0x3;
	unsigned int c5 = (kmer >> 8) & 0x3;

	std::ostringstream file_name;
	file_name << bitstoc(c1) << bitstoc(c2) << bitstoc(c3) << bitstoc(c4) << bitstoc(c5);
	return file_name.str();
}

std::vector<KmerFile *> FastqReader::generateKmers(size_t kmersize)
{
	std::map<std::string, std::vector<uint64_t>> buffer;
	size_t counter = 0;
	std::string line;
	while (readNextSequence(line))
	{
		auto kmers = findKmers(line, kmersize);
		for (uint64_t kmer : kmers)
		{
			auto firstFive = getfirstFiveChars(kmer);
#ifdef USE_GZIP 
			buffer[firstFive].push_back(kmer);
			counter++;
			if (counter == 1000)
			{
				for (auto &pair : buffer)
				{
					if (pair.second.size() > 0)
					{
						writeCompressedFile(filename_map[pair.first], &pair.second[0], pair.second.size() * sizeof(uint64_t));
					}
				}
				buffer.clear();
			}
#else
			*filename_map[firstFive] << kmer << '\n';
#endif
			if (filename_fastq_map[firstFive])
			{
				filename_fastq_map[firstFive]->kmercount++;
			} 
			else
			{
				auto fastqFile = new KmerFile();
				fastqFile->filename = std::string("temp/" + firstFive + ".txt").c_str();
				fastqFile->kmercount = 1;
				filename_fastq_map[firstFive] = fastqFile;
			}
		}
#ifdef USE_GZIP 
			for (auto &pair : buffer)
			{
				if (pair.second.size() > 0)
				{
					writeCompressedFile(filename_map[pair.first], &pair.second[0], pair.second.size() * sizeof(uint64_t));
				}
			}
			buffer.clear();
#endif
		/*
		 *
		auto sequences = findKmers(line, kmersize);
		for (auto &sequence : sequences)
		{
			// FIXME replace N
			//sequence.replace("N", "A");
			auto seq = sequence.substr(0, 5);
			if (sequence.find('N') == std::string::npos)
			{
				buffer_map[seq] << sequence << '\n';
			}
		}
		buffer_counter++;
		if (buffer_counter == 20)
		{
			insert_sequence(buffer_map);
			buffer_map.clear();
			buffer_counter = 0;
		}
		 */
	}

	std::vector<KmerFile *> files;
	for (auto &pair : filename_fastq_map)
	{
		files.push_back(pair.second);
	}
	for (auto &pair : filename_map)
	{
#ifdef USE_GZIP 
		closeGzipFile(pair.second);
#else
		if (pair.second->is_open())
		{
			pair.second->close();
			delete pair.second;
		}
#endif
	}
	filename_fastq_map.clear();
	filename_map.clear();

	if (m_file.is_open())
	{
		m_file.close();
	}
	return files;
}

bool FastqReader::readNextSequence(std::string &line)
{
	/*
	 * Assumed that, the structure of file is OK
	 * Therefore no defensive code written
	 */
	std::string temp;
	m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	//std::getline(*m_file, temp); // skip sequence id line 
	auto b = std::getline(m_file, line) ? true : false; 
	//std::getline(*m_file, temp); // skip + line
	m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	//std::getline(*m_file, temp);  // skip quality line
	return b;
}
