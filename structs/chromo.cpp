/*
 *
 */
#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>

#include "chromo.h"


Chromo::Chromo (std::string _filename, std::string _description, size_t _start, size_t _end): 
  filename(_filename), 
  description(_description),
  start(_start),
  end(_end),
  sequence(NULL)
{
  for (size_t i=0; i < (end - start)/READ_BLOCK_SIZE + 1; i++) {
    //NOTE: it's quite dirty way to avoid sign-error for size_t.
    size_t start_offset;
    if (start + READ_BLOCK_SIZE * i < INTERSECTION_SIZE) {
      start_offset = 0;
    } else {
      start_offset = start + READ_BLOCK_SIZE * i - INTERSECTION_SIZE;
    }
    
    start_part.push_back(std::max(start, start_offset ));
    end_part.push_back(std::min(end, start + READ_BLOCK_SIZE*(i+1)));
    length_part.push_back(end_part[i] - start_part[i]);
  }
  for (size_t i=0; i < (end - start)/READ_BLOCK_SIZE + 1; i++) {
    std::cout << "Parts: " <<start_part[i] << "   " << end_part[i] << "    " << length_part[i] << std::endl;
  }

}

std::string& Chromo::getFilename() {
  return filename;
}

std::string& Chromo::getDescription() {
  return description;
}


char* Chromo::getSeqPtr(size_t part_number) {
  /*    sequence = new char[end-start+1];
    std::ifstream fin(filename.c_str(), std::ifstream::in);
    fin.seekg(start);
    fin.read(sequence, end-start);
    sequence[end-start] = '\0';
    length = end - start; */
    if (part_number > this->getNumberOfParts() - 1) {
      throw std::invalid_argument("Incorrect part_number");
    }
  
    if (sequence != NULL) {
      delete[] sequence;
      sequence = NULL;
    }
    sequence = new char[end_part[part_number] - start_part[part_number] + 1];
    std::ifstream fin(filename.c_str(), std::ifstream::in);
    fin.seekg(start_part[part_number]);
    fin.read(sequence, end_part[part_number] - start_part[part_number]);
    sequence[end_part[part_number] - start_part[part_number]] = '\0';
    
    return sequence;
}

//This returns absolute offset of given part in file.
size_t Chromo::getPartOffset ( size_t part_number ) {
    if (part_number > this->getNumberOfParts() - 1) {
      throw std::invalid_argument("Incorrect part_number");
    }
    return start_part[part_number];
}

//TODO: std::vector<char> getWord (position, length);

size_t Chromo::getNumberOfParts() {
    return (end - start)/READ_BLOCK_SIZE + 1;
}

size_t Chromo::getPartLength(size_t part_number) {
    return length_part[part_number];
}

void Chromo::releasePart ( size_t part_number ) {
    if (sequence != NULL) {
      delete[] sequence;
      sequence = NULL;
    }
}




Chromo::~Chromo () {
  if (sequence != NULL) {
    delete[] sequence;
  }
}
