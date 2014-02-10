#ifndef CHROMO_H
#define CHROMO_H

#include <string>


class Chromo {
  std::string filename;
  std::string description;
  std::size_t start;
  std::size_t end;
  char * sequence;
  std::size_t length;
  
public:
    Chromo (std::string _filename, std::string _description, std::size_t _start, std::size_t _end);
    
    std::string& getFilename();
    std::string& getDescription();
    
    char* getSeqPtr();
// NOTE: We are returning size of chromo as if it was stored in string: without null character
    size_t size();
   ~Chromo ();
};

#endif // CHROMO_H
