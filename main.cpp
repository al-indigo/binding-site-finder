#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <stack>
#include <time.h>

//#include <boost/regex.hpp>
#include "ahoc/AhoCorasickPlus.h"
#include "structs/pwm.h"

#define DTIME_PROFILING

int main(int argc, char **argv) {

#ifdef DTIME_PROFILING
  double tstart, tstop, ttime;
  tstart = (double)clock()/CLOCKS_PER_SEC;
#endif

  std::ifstream ifs("/Volumes/Data/Downloads/genome/chr1.plain");
  std::string seq( (std::istreambuf_iterator<char>(ifs) ), (std::istreambuf_iterator<char>() ) );

#ifdef DTIME_PROFILING
  tstop = (double)clock()/CLOCKS_PER_SEC;
  std::cout << "File has been read in " << tstop - tstart << " seconds" << std::endl; 
#endif

#ifdef DTIME_PROFILING
  tstart = (double)clock()/CLOCKS_PER_SEC;
#endif

  std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
#ifdef DTIME_PROFILING
  tstop = (double)clock()/CLOCKS_PER_SEC;
  std::cout << "File has been upper-cased in " << tstop - tstart << " seconds" << std::endl;
#endif

#ifdef DTIME_PROFILING
  tstart = (double)clock()/CLOCKS_PER_SEC;
#endif

//  Pwm matrix("/Users/al/Downloads/pwms/AHR_si.pat", 8);
  Pwm matrix("/Users/al/Downloads/pwms/TAL1_f2.pat", 16.6);

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
    patterns = matrix.getWords(16.6, 50000);
    for (size_t i = 0; i < patterns.size(); i++)
    {
        AhoCorasickPlus::EnumReturnStatus status;
        AhoCorasickPlus::PatternId patId = totalWords++;
  //      status = atm.addPattern(patterns.top(), patId);
  //      patterns.pop();
        status = atm.addPattern(patterns[i], patId);
        if (status!=AhoCorasickPlus::RETURNSTATUS_SUCCESS)
        {
            std::cout << "Failed to add: " << std::endl;
        }
    }
  }
  atm.finalize();

#ifdef DTIME_PROFILING
  tstop = (double)clock()/CLOCKS_PER_SEC;
  std::cout << "Automata is initialized in " << tstop - tstart << " seconds"  << std::endl;
#endif

#ifdef DTIME_PROFILING
  tstart = (double)clock()/CLOCKS_PER_SEC;
#endif

  AhoCorasickPlus::Match aMatch;
  size_t occurances = 0;
  atm.search(seq, true);
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

  return 0;
}
