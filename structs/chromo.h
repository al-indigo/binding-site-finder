#ifndef CHROMO_H
#define CHROMO_H

/*NOTE: Important! I assume here that there is no words longer than 128 symbols.
 *      In fact I'm searching the same places but we getting rid of duplicates in merge-sort.
 */
#define INTERSECTION_SIZE 128LLU

/*NOTE: Here the part-size that is stored in-memory is defined. 
 *      You may change it if you want but I've measured -- it's optimal 
 *      for regular HDD.*/
#define READ_BLOCK_SIZE size_t(32*1024*1024LLU - INTERSECTION_SIZE)


#include <string>
#include <vector>

class Chromo {
  std::string filename;
  std::string description;
  std::size_t start;
  std::size_t end;
  std::vector<std::size_t> start_part;
  std::vector<std::size_t> end_part;
  std::vector<std::size_t> length_part;
  char * sequence;
  std::size_t length;
  
public:
    Chromo (std::string _filename, std::string _description, std::size_t _start, std::size_t _end);
    
    std::string& getFilename();
    std::string& getDescription();
    size_t getNumberOfParts();
    
//NOTE: part number start from 0
    char* getSeqPtr(size_t part_number);
    size_t getPartOffset(size_t part_number);
    
// NOTE: We are returning size of chromo as if it was stored in string: without null character
    size_t getPartLength (size_t part_number);
    void releasePart (size_t part_number);
   ~Chromo ();
};

#endif // CHROMO_H
