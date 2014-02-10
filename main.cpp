#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <stack>
#include <time.h>

//#include <boost/regex.hpp>
#include "ahoc/AhoCorasickPlus.h"
#include "structs/pwm.h"
#include "structs/chromovector.h"

#define DTIME_PROFILING

int main(int argc, char **argv) {
  
  size_t test_start = 2400;
  size_t test_end = 50000;

  std::vector<std::string> filenames;
  std::vector<std::string> descriptions;
  std::vector<size_t> starts;
  std::vector<size_t> ends;
  
  filenames.push_back(std::string("/Users/al/Downloads/genome/chr1_upper.plain"));
  descriptions.push_back(std::string("chr1"));
  starts.push_back(test_start);
  ends.push_back(test_end);
  
  filenames.push_back(std::string("/Users/al/Downloads/genome/chr1_upper.plain"));
  descriptions.push_back(std::string("chr1"));
  starts.push_back(test_start);
  ends.push_back(test_end);
  
  ChromoVector sequences(filenames, descriptions, starts, ends);

#ifdef DTIME_PROFILING
  double tstart, tstop, ttime;
#endif

  /*
  try {} catch (std::exception e) {
#ifdef DDEBUG_PRINT
    std::cerr << e.what();
#endif
  }
*/
  
#ifdef DTIME_PROFILING
  tstart = (double)clock()/CLOCKS_PER_SEC;
#endif

  Pwm matrix("/Users/al/Downloads/pwms/AHR_si.pat", 3);
//  Pwm matrix("/Users/al/Downloads/pwms/TAL1_f2.pat", 13.6);

//  std::stack<std::string> patterns = matrix.getWords(8, 10000000);;
//  while(!matrix.hasMoreWords()) patterns = matrix.getWords(8, 10000000);
  STRING_CONTAINER patterns;
//  patterns.push_back("AGATAATCTTTATTGTCC");
//  patterns = matrix.getWords(8, 10000000);

  AhoCorasickPlus atm;
  
//  std::cout << "#patterns fit: " << patterns.size() << std::endl;
  size_t totalWords = 0;
  while (!matrix.hasMoreWords()) {
    STRING_CONTAINER patterns;
    patterns = matrix.getWords(3, 50000);
    for (size_t i = 0; i < patterns.size(); i++)
    {
        AhoCorasickPlus::EnumReturnStatus status;
        AhoCorasickPlus::PatternId patId = totalWords++;
  //      status = atm.addPattern(patterns.top(), patId);
  //      patterns.pop();
        status = atm.addPattern(patterns[i], patId);
        if (status!=AhoCorasickPlus::RETURNSTATUS_SUCCESS)
        {
            std::cerr << "Failed to add: " << std::endl;
        }
    }
  }
  atm.finalize();
  
  for (size_t s_i = 0; s_i < sequences.size(); s_i++) {
    char * seq = sequences.getSeq(s_i);
    size_t length = sequences.getSeqLength(s_i);

  #ifdef DTIME_PROFILING
    tstop = (double)clock()/CLOCKS_PER_SEC;
    std::cout << "Automata is initialized in " << tstop - tstart << " seconds"  << std::endl;
  #endif

  #ifdef DTIME_PROFILING
    tstart = (double)clock()/CLOCKS_PER_SEC;
  #endif

    AhoCorasickPlus::Match aMatch;
    size_t occurances = 0;
    std::cout << "Results for " << sequences.getFilename(s_i) << std::endl;
    atm.search(seq, length, false);
    while (atm.findNext(aMatch))
    {
  //    std::cout << "@" << aMatch.position << "\t#" << aMatch.id << "\t" << sample_patterns[aMatch.id] << std::endl;
      std::cout << aMatch.position << std::endl;
      occurances++;
    }
    
  #ifdef DTIME_PROFILING
    tstop = (double)clock()/CLOCKS_PER_SEC;  
    std::cout << "Search performed in " << tstop - tstart << " seconds" << std::endl;
  #endif
    std::cout << "Occurances found: " << occurances << std::endl;
  }

  return 0;
}
