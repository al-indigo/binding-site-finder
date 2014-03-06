#ifndef MERGESORT_INFILE_H
#define MERGESORT_INFILE_H

#include <string>

#define FSTREAM_BUF_SIZE 16*1024*1024LLU

void merge_sort(std::string filename1in, std::string filename2in, std::string filename_out);

#endif