/*
 *
 */
#include <exception>
#include <iostream>
#include <fstream>

#include "chromo.h"


Chromo::Chromo (std::string _filename, std::string _description, size_t _start, size_t _end): 
  filename(_filename), 
  description(_description),
  start(_start),
  end(_end)
{
    sequence = new char[end-start+1];
    std::ifstream fin(filename.c_str(), std::ifstream::in);
    fin.seekg(start);
    fin.read(sequence, end-start);
    sequence[end-start] = '\0';
    length = end - start;
}

std::string& Chromo::getFilename() {
  return filename;
}

std::string& Chromo::getDescription() {
  return description;
}


char* Chromo::getSeqPtr() {
    return sequence;
}

size_t Chromo::size() {
    return length;
}


Chromo::~Chromo () {
  
  delete[] sequence;

}
