#ifndef CHROMO_H
#define CHROMO_H

/* IMPORTANT: more block you use, faster the program works -- it needs to make less
 *            merges of results which take lot of time. Block size is in bytes
 */
#define READ_BLOCK_SIZE size_t(96*1024*1024LLU)


#include <string>
#include <vector>
#include <set>
#include <mutex>
#include <atomic>
#include "pwm.h"

class Chromo {
  std::string filename;
  std::string description;
  std::size_t start;
  std::size_t end;
  
  size_t intersection_size; //NOTE: intersection size must be equal to length of matrix, otherwise in multithreaded mode there will be race-condition.
//  size_t block_size;        //It's real read block: it's equal READ_BLOCK_SIZE-intersection_size
  
  std::vector<std::size_t> start_part;
  std::vector<std::size_t> end_part;
  std::vector<std::size_t> length_part;
  char * sequence;
  size_t seq_start;
  std::atomic_size_t current_part;
  size_t number_of_parts;

  std::mutex lock;
#ifdef DDEBUG_PRINT
  std::mutex print_lock;
#endif
  
public:
    Chromo (std::string _filename, std::string _description, std::size_t _start, std::size_t _end, size_t _intersection_size);
    
    std::string& getFilename();
    std::string& getDescription();
    size_t getNumberOfParts();
    
//NOTE: part number start from 0
    char* getSeqPtr(size_t part_number);
    size_t getPartOffset(size_t part_number);
    
// NOTE: We are returning size of chromo as if it was stored in string: without null character
    size_t getPartLength (size_t part_number);
    void releasePart (size_t part_number);
    
    void getWordsAsPaths (std::vector<size_t>& positions, size_t length, std::vector <std::vector<char> >& result);
    void getWordAsPathTest (size_t position, size_t length, std::vector<char>& result);
    bool getWordScores (size_t position, Pwm& matrix, std::pair<double, double>& scores);
    
   ~Chromo ();
};

#endif // CHROMO_H
