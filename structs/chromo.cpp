/*
 *
 */
#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <assert.h>
#include <cstring>
#include <cstdio>
#include "chromo.h"


Chromo::Chromo (std::string _filename, std::string _description, size_t _start, size_t _end): 
  filename(_filename), 
  description(_description),
  start(_start),
  end(_end),
  sequence(NULL),
  current_part(-1)
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
    if (part_number > this->getNumberOfParts() - 1) {
      throw std::invalid_argument("Incorrect part_number");
    }
  
    if (sequence != NULL) {
      delete[] sequence;
      sequence = NULL;
    }
    sequence = new char[end_part[part_number] - start_part[part_number] + 1];
//    std::ifstream fin(filename.c_str(), std::ifstream::in);
    FILE * fin = fopen(filename.c_str(), "r");
//    fin.seekg(start_part[part_number]);
    fseek(fin, start_part[part_number], SEEK_SET);
//    fin.read(sequence, end_part[part_number] - start_part[part_number]);
    fread(sequence, sizeof(char), end_part[part_number] - start_part[part_number], fin);
    sequence[end_part[part_number] - start_part[part_number]] = '\0';
    current_part = part_number;
    
    return sequence;
}



//This returns absolute offset of given part in file.
size_t Chromo::getPartOffset ( size_t part_number ) {
    if (part_number > this->getNumberOfParts() - 1) {
      throw std::invalid_argument("Incorrect part_number");
    }
    return start_part[part_number];
}

inline static void wordToPath (char * buffer, char * destination, size_t length) {
  for (size_t i = 0; i < length; i++) {
    if (buffer[i] == 'A') {
        destination[i] = 0;
    } else if (buffer[i] == 'C') {
        destination[i] = 1;
    } else if (buffer[i] == 'G') {
        destination[i] = 2;
    } else if (buffer[i] == 'T') {
        destination[i] = 3;
    } else {
      assert(true);
    }
  }
}


void Chromo::getWordsAsPaths (std::vector<size_t>& positions, size_t length, std::vector <std::vector<char> >& result) {
//  for (std::vector<std::vector<char> >::iterator i = result.begin(); i != result.end(); ++i) {
//    (*i).resize(length, -1);
//  }
  
//  std::cout << "Preparing to put words, size reserved is " << result[0].size() << " model len is " << length << " need to put " << positions.size() << " words as byte paths" << std::endl;
  
  std::ifstream fin(filename.c_str(), std::ifstream::in);

  char buffer[length];
  size_t counter = 0;
  for (std::vector<size_t>::iterator i = positions.begin(); i != positions.end(); ++i, counter++) {
    fin.seekg(*i);
    fin.read(buffer, length);
    wordToPath(buffer, &((result[counter])[0]), length);
  }
  fin.close();
}

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
    current_part = -1;
}

void Chromo::getWordAsPathTest (size_t position, size_t length, std::vector<char>& result) {
  //guess part, where word is located
  size_t part_no = 0;
  if (position >= end_part[getNumberOfParts()-1] - length) {
    std::cerr << "This should never happen! Bounds violation" << std::endl;
    exit(-1);
  }
  
  result.resize(length);
  
  for (size_t i = 0; i < end_part.size(); i++) {
    if (position < end_part[i] - length) {
      part_no = i;
      break;
    }
  }
  if (current_part != part_no) {
    sequence = getSeqPtr(part_no);
  }
  
  size_t read_start = position - getPartOffset(current_part);
//  std::cout << read_start << "\t"  << length << std::endl;
  char buffer[length];
  memcpy(buffer, sequence + read_start, length);
  wordToPath(buffer, &(result[0]), length);
 
}

// size_t Chromo::getFullSequenceAsPaths (size_t stop_point, size_t length, std::vector <std::vector<char> >& result) {
//   sequence = getSeqPtr(part_number); // I know we may re-read this part accidentally, but it's much more stable way.
//   size_t cur_len = getPartLength(part_number);
//   
//   
//   
//   
// }




Chromo::~Chromo () {
  if (sequence != NULL) {
    delete[] sequence;
  }
}
