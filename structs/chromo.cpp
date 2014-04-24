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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "chromo.h"

Chromo::Chromo (std::string _filename, std::string _description, size_t _start, size_t _end, size_t _intersection_size): 
  filename(_filename), 
  description(_description),
  start(_start),
  end(_end),
  intersection_size(_intersection_size),
  sequence(NULL),
  seq_start(0),
  current_part(-1)
{
  char * buf = new char[READ_BLOCK_SIZE];
  FILE * fin = fopen(filename.c_str(), "r");
  fseek(fin, 0, SEEK_SET);
  size_t bytes_read = fread(buf, sizeof(char), READ_BLOCK_SIZE, fin);
  if (bytes_read == 0) {
    std::cerr << "Chromosome file " << filename << " is empty, can not continue" << std::endl;
    exit(-1);
  }
  // If it's fasta with description, we'll parse it here.
  if (buf[0] == '>') {
    char * linebr = strchr(buf, '\n');
    if (linebr == NULL) {
      seq_start = 0;
    } else {
      seq_start = linebr - buf + 1;
    }
  }
  fclose(fin);
  
  size_t full_sequence_length = 0;
  struct stat filestatus;
  stat(filename.c_str(), &filestatus);
  full_sequence_length = filestatus.st_size - seq_start;
  
  end = std::min(end, full_sequence_length);
  if (end == 0) { end = full_sequence_length; }
  delete[] buf;
  
//  std::cout << "end: " << end << " file size: " << filestatus.st_size << " full seq: " << full_sequence_length << " offset: " << seq_start <<"\n" ;
  
  for (size_t i=0; i < (end - start)/READ_BLOCK_SIZE + 1; i++) {
    size_t start_offset;
    if (start + READ_BLOCK_SIZE * i < intersection_size) {
      start_offset = 0;
    } else {
      start_offset = start + READ_BLOCK_SIZE * i - intersection_size;
    }
    
    start_part.push_back(std::max(start, start_offset ));
    end_part.push_back(std::min(end, start + READ_BLOCK_SIZE*(i+1)));
    length_part.push_back(end_part[i] - start_part[i]);
  }
#ifdef DDEBUG_PRINT  
  for (size_t i=0; i < (end - start)/READ_BLOCK_SIZE + 1; i++) {
    std::cout << "Parts: " <<start_part[i] << "   " << end_part[i] << "    " << length_part[i] << std::endl;
  }
#endif
  number_of_parts = (end - start)/READ_BLOCK_SIZE + 1;
}

std::string& Chromo::getFilename() {
  return filename;
}

std::string& Chromo::getDescription() {
  return description;
}


char* Chromo::getSeqPtr(size_t part_number) {
    
    
    if (part_number > number_of_parts - 1) {
      throw std::invalid_argument("Incorrect part_number");
    }
    
    lock.lock();
#ifdef DDEBUG_PRINT    
    std::cout << "Locked by:" << lock.native_handle() << ", current_part: " << current_part << ", part to acquire:" << part_number << " I'm thread " << std::this_thread::get_id() << std::endl;
#endif
    if (current_part.load() != part_number) {
      releasePart(current_part);
    } else {
      lock.unlock();
      return sequence;
    }
    
    size_t plen = end_part[part_number] - start_part[part_number];
    
    sequence = new char[plen + 1 + SUPERALPHABET_SIZE];
    memset(sequence, 0, plen + 1 + SUPERALPHABET_SIZE);
    FILE * fin = fopen(filename.c_str(), "r");
    fseek(fin, seq_start + start_part[part_number], SEEK_SET);
    size_t bytes_read = fread(sequence, sizeof(char), plen, fin);
    if (bytes_read == 0) {
      std::cerr << "Read error occured while reading chromosome file. Can not continue." << std::endl;
      exit(-1);
    }
    
    for (auto i = 0; i < plen; i++) {
      sequence[i] &=0xDF; //uppercase trick
    }
    
//    sequence[plen] = '\0';
    current_part.store(part_number);
    
    fclose(fin);
    lock.unlock();
    return sequence;
}



//This returns absolute offset of given part in file.
size_t Chromo::getPartOffset ( size_t part_number ) {
    if (part_number > this->number_of_parts - 1) {
      throw std::invalid_argument("Incorrect part_number");
    }
    return start_part[part_number];
}


inline static bool wordToPath (char * buffer, char * destination, size_t length) {
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
      return false; //NOTE: in most cases that means that we are at "N" letter.
    }
  }
  return true;
}

inline static bool wordToPathBit (char * buffer, char * destination, size_t length) {
  bool res = 0x00;
  
  for (size_t i=0; i< length; i++) {
    char bitrepr = 0x0F;
    bitrepr &= buffer[i];
    res |= (buffer[i]&0x08)>>3;
    destination[i] = (bitrepr&0x07)>>1;
  }
  return !res;
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
    return number_of_parts;
}

size_t Chromo::getPartLength(size_t part_number) {
    return length_part[part_number];
}

void Chromo::releasePart ( size_t part_number ) {
  if (sequence != NULL) {
    delete[] sequence;
    sequence = NULL;
  }
  current_part.store(-1);
}

void Chromo::getWordAsPathTest (size_t position, size_t length, std::vector<char>& result) {
  //guess part, where word is located
  size_t part_no = 0;
  if (position >= end_part[number_of_parts-1] - length) {
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
  if (current_part.load() != part_no) {
    sequence = getSeqPtr(part_no);
  }
  
  size_t read_start = position - getPartOffset(current_part.load());
//  std::cout << read_start << "\t"  << length << std::endl;
  char buffer[length];
  memcpy(buffer, sequence + read_start, length);
  wordToPath(buffer, &(result[0]), length);
 
}

bool Chromo::getWordScores (size_t position, Pwm& matrix, std::pair<double, double>& scores) {
  //guess part, where word is located
  size_t part_no = 0;
  if (position >= end_part[number_of_parts-1] - intersection_size) {
    std::cerr << "This should never happen! Bounds violation" << std::endl;
    exit(-1);
  }
  
  for (size_t i = 0; i < end_part.size(); i++) {
    if (position < end_part[i] - intersection_size) {
      part_no = i;
      break;
    }
  }
  
  if (current_part.load() != part_no) {
#ifdef DDEBUG_PRINT    
    print_lock.lock();
    std::cout << "Want to change part from " << current_part.load() << " to " << part_no << ", I'm at " << position << ", I'm thread " << std::this_thread::get_id() << std::endl;
    print_lock.unlock();
#endif
    sequence = getSeqPtr(part_no);
  }

  size_t read_start = position - start_part[current_part.load()];
//  std::cout << read_start << "\t"  << length << std::endl;
  char buffer[intersection_size+SUPERALPHABET_SIZE];
  memset(buffer, 0, intersection_size+SUPERALPHABET_SIZE);
  if (matrix.getMatrixType() == pat) {
    if(!wordToPath(sequence + read_start, buffer, intersection_size)) {
      return false;
    }
  } else if (matrix.getMatrixType() == pat_bit_optimization) {
    if(!wordToPathBit(sequence + read_start, buffer, intersection_size)) {
      return false;
    }
  }
  
  return matrix.getScorePair(buffer, scores);
}


void Chromo::getManyWordScores (std::vector<size_t>& positions_to_read_again, Pwm& matrix, std::vector<double>& scores, std::vector<size_t>& matched, std::vector<bool>& strand) {
  if (positions_to_read_again[positions_to_read_again.size()-1] >= end_part[number_of_parts-1] - intersection_size) {
    std::cerr << "This should never happen! Bounds violation" << std::endl;
    exit(-1);
  }
  char buffer[intersection_size+SUPERALPHABET_SIZE];
  for (auto&& position: positions_to_read_again) {
    size_t part_no = 0;
    
    for (size_t i = 0; i < end_part.size(); i++) {
      if (position < end_part[i] - intersection_size) {
        part_no = i;
        break;
      }
    }
    
    if (current_part.load() != part_no) {
  #ifdef DDEBUG_PRINT    
      print_lock.lock();
      std::cout << "Want to change part from " << current_part.load() << " to " << part_no << ", I'm at " << position << ", I'm thread " << std::this_thread::get_id() << std::endl;
      print_lock.unlock();
  #endif
      sequence = getSeqPtr(part_no);
    }

    size_t read_start = position - start_part[current_part.load()];

    memset(buffer, 0, intersection_size+SUPERALPHABET_SIZE);
    if(!wordToPathBit(sequence + read_start, buffer, intersection_size)) {
      continue;
    }
    matrix.getScores(buffer, scores, strand, matched, position);
  }
}


std::vector<bool> wordToPathBitVector (char * buffer, char * destination, size_t length, size_t word_length) {
  std::vector<bool> is_invalid_marker(length+1, true);
  long int wlen = (int) word_length;
  std::vector<bool> need_to_check_bool(length+1, true);
  std::vector<int> need_to_ignore_neighbourhood;
  
  for (size_t i=0; i< length; i++) {
    char bitrepr = 0x0F;
    bitrepr &= buffer[i];
    is_invalid_marker[i] = (0x00|(buffer[i]&0x08)>>3);
    destination[i] = (bitrepr&0x07)>>1;
  }

  for (int i = is_invalid_marker.size() - 1; i >= 0; i--) {
    if (is_invalid_marker[i]) {
      for (int k = 0; k <= wlen && i > 0; k++) {
        need_to_check_bool[i] = false;
        i--;
        if (is_invalid_marker[i]) {
          need_to_check_bool[i] = false;
          k = 0;
        }
      }
    }
  }

  return need_to_check_bool;
}


void Chromo::getWordScoresVector (size_t first, size_t last, Pwm& matrix, std::vector<double>& scores, std::vector<bool>& strand, std::vector<size_t>& positions) {
  //guess part, where word is located
  size_t part_no = 0;
  if (last > end_part[number_of_parts-1]) {
    std::cerr << "This should never happen! Bounds violation: last " << last << ", end: " << end_part[number_of_parts-1] - intersection_size << std::endl;
    exit(-1);
  }
  
  for (size_t i = 0; i < end_part.size(); i++) {
    if (first < end_part[i] - intersection_size) {
      part_no = i;
      break;
    }
  }
  
  if (current_part.load() != part_no) {
    sequence = getSeqPtr(part_no);
  }

  size_t read_start = first - start_part[current_part.load()];

  char * buffer = new char[last - first + intersection_size];

  
  std::vector<bool> need_to_check = wordToPathBitVector(sequence + read_start, buffer, last-first+intersection_size, intersection_size);
  matrix.getScoresVector(buffer, need_to_check, scores, strand, positions, first);

  delete [] buffer;
  return;
}


Chromo::~Chromo () {
  if (sequence != NULL) {
    delete[] sequence;
  }
}
