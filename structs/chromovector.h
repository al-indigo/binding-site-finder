#ifndef CHROMOVECTOR_H
#define CHROMOVECTOR_H

#include "chromo.h"

#include <vector>

class ChromoVector {
  std::vector<Chromo*> chromovector;
  
public:
  ChromoVector (std::vector<std::string> filenames, 
                std::vector<std::string> description, 
                std::vector<size_t> start, 
                std::vector<size_t> end);
  size_t size();
  std::string &getFilename(size_t sequence_number);
  std::string &getDescription(size_t sequence_number);
  char * getSeq (size_t sequence_number);
  size_t getSeqLength (size_t sequence_number);
  
 ~ChromoVector ( );
};

#endif // CHROMOVECTOR_H
