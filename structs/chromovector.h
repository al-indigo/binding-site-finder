#ifndef CHROMOVECTOR_H
#define CHROMOVECTOR_H

#include "chromo.h"

#include <vector>

class ChromoVector {
  std::vector<Chromo*> chromovector;
  std::vector<std::vector<std::string> > tempresults;
  std::vector<std::string> result_filenames;

public:
  ChromoVector (std::vector<std::string> filenames,
                std::vector<std::string> description,
                std::vector<size_t> start,
                std::vector<size_t> end);
  size_t size();
  std::string& getFilename(size_t sequence_number);
  std::string& getDescription(size_t sequence_number);
  char * getSeq (size_t sequence_number, size_t part_number);
  size_t getNumberOfParts (size_t sequence_number);
  size_t getPartLength ( size_t sequence_number, size_t part_number );
  size_t getSeqLength (size_t sequence_number);
  size_t getAbsoluteOffset (size_t sequence_number, size_t part_number);
  void releasePart ( size_t sequence_number, size_t part_number);
  
  void getWordsAsPaths (size_t sequence_number, std::vector<size_t>& positions, size_t length, std::vector <std::vector<char> >& result) { chromovector[sequence_number]->getWordsAsPaths(positions, length, result); };
  void getWordAsPathTest (size_t sequence_number, size_t position, size_t length, std::vector<char>& result) { chromovector[sequence_number]->getWordAsPathTest(position, length, result); }
  
 ~ChromoVector ( );
};

#endif // CHROMOVECTOR_H
