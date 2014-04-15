/*
 *
 */

#include "chromovector.h"
#include <stdexcept>

ChromoVector::ChromoVector (std::vector<std::string> filenames, 
                            std::vector<std::string> description,
                            std::vector<size_t> start, 
                            std::vector<size_t> end,
                            size_t intersection_size
                           ) {
  if (!(filenames.size() == description.size() && 
        description.size() == start.size() &&
        start.size() == end.size() &&
        end.size() > 0)) {
    throw std::invalid_argument("Incorrect number of arguments");
  } else {
    chromovector.reserve(filenames.size());
    for (size_t i = 0; i < filenames.size(); i++) {
      chromovector.push_back(new Chromo(filenames[i], description[i], start[i], end[i], intersection_size));
    }
  }
}


size_t ChromoVector::size() {
  return chromovector.size();
}

std::string& ChromoVector::getFilename ( size_t sequence_number ) {
  return chromovector[sequence_number]->getFilename();
}

std::string& ChromoVector::getDescription ( size_t sequence_number ) {
  return chromovector[sequence_number]->getDescription();
}

char* ChromoVector::getSeq ( size_t sequence_number, size_t part_number ) {
  return chromovector[sequence_number]->getSeqPtr(part_number);
}

size_t ChromoVector::getPartLength ( size_t sequence_number, size_t part_number ) {
  return chromovector[sequence_number]->getPartLength(part_number);
}

size_t ChromoVector::getNumberOfParts ( size_t sequence_number ) {
  return chromovector[sequence_number]->getNumberOfParts();
}

void ChromoVector::releasePart ( size_t sequence_number, size_t part_number ) {
  chromovector[sequence_number]->releasePart(part_number);
}

size_t ChromoVector::getAbsoluteOffset ( size_t sequence_number, size_t part_number ) {
  return chromovector[sequence_number]->getPartOffset(part_number);
}



ChromoVector::~ChromoVector () {
  for (size_t i = 0; i < this->size(); i++) {
    delete chromovector[i];
  }
}
