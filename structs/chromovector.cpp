/*
 *
 */

#include "chromovector.h"

ChromoVector::ChromoVector (std::vector<std::string> filenames, 
                            std::vector<std::string> description, 
                            std::vector<size_t> start, 
                            std::vector<size_t> end) {
  if (!(filenames.size() == description.size() && 
        description.size() == start.size() &&
        start.size() == end.size() &&
        end.size() > 0)) {
    throw std::invalid_argument("Incorrect number of arguments");
  } else {
    for (size_t i = 0; i < filenames.size(); i++) {
      chromovector.push_back(new Chromo(filenames[i], description[i], start[i], end[i]));
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


char* ChromoVector::getSeq ( size_t sequence_number ) {
  return chromovector[sequence_number]->getSeqPtr();

}

size_t ChromoVector::getSeqLength ( size_t sequence_number ) {
  return chromovector[sequence_number]->size();
}


ChromoVector::~ChromoVector () {
  for (size_t i = 0; i < this->size(); i++) {
    delete chromovector[i];
  }
}
