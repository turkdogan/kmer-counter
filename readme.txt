INTRODUCTION
* A Sequence diagram is generated for quick reference (sequence_diagram.png)
* The k-mer counter is written in C++ 11. 
* It is tested on Windows and Linux.
* I tried to use C++ syntax/functions over C functions, which may result a bit slow execution.

RUNNING
* To compile and run:
 - make
 - ./kmer_counter --filename example.fastq --kmersize 30 --topcount 25
 - Optional --outfilename {filename} can be given, or else out.txt is used

ZLIB
* Zlib library is enabled and required for Linux, but FASTQ file MUST BE "uncompressed" before given as parameter.
* Zlib is not required for windows, 10x fastq-size disk is required.
 - Even compressed, at least 3x of fastq file size is required for linux.
* Zlib is enabled default for Linux, this increases running time, decreases disk time
* To enable/disable gzip uncomment/comment USE_GZIP definition line in the common.h file

TECHNICAL DETAILS
* The Fastq file is preprocessed and splitted into buckets,
* Each bucket contains specific kmers categorized by the first three letters of the kmer
* Each bucket processing is carried out isolated from other buckets,
 - Therefore each bucket processing is performed in a seperate thread,
 - As a result, the program can be extended to a parallel/distributed structure easily.
* I focused on memory restriction, as a result; performance is not academic grade.
* A Bloom filter is used to process each bucket.
* Bloom filter uses Murmur hashing algorithm.
* There is room to improve Bloom filter implementation, especially the bitset structure.
* Murmur hashing implementation is grabbed from Internet,
 - It is one of the most suitable hashing algorithm for Bloom Filter.
* Codes are tested 64 bit host computer. I used 64 bit version of Murmur hashing algorithm
* The sequences containing 'N' are ignored. 
 - Due to lack of domain info, I decided to ignore istead of replacing 'N' with '{ACGT}'

UNIT TESTS
* Tests folder includes unit tests
* Catch library is used for the unit tests

FINAL NOTES/KNOWN ISSUES
* Not have enough time for the following issues
 - Open FastQ.gz directly without
 - Processing sequences containing 'N'
 - Mapping sequence to long such that A => 00, C = 01, G = 10, T = 11
 - MacOS is not tested
