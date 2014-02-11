#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <stack>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "ahoc/AhoCorasickPlus.h"
#include "structs/pwm.h"
#include "structs/chromovector.h"
#include "utils/mergesort_infile.h"

#define DTIME_PROFILING

int main(int argc, char **argv) {
//NOTE: results_path should exist (I'm ruling all the filenames stuff with python wrapper - it's easier)
  std::string result_folder("/Users/al/Programming/binding-site-finder/build/results/");
 
  size_t test_start = 0;
  size_t test_end = 100;
  
  double score_threshold = -5.0;
  
  std::vector<std::string> filenames;
  std::vector<std::string> descriptions;
  std::vector<std::string> result_filenames;
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
  
  result_filenames.push_back(std::string("chr1-result"));
  result_filenames.push_back(std::string("chr2-result"));
  
  ChromoVector sequences(filenames, descriptions, result_filenames, starts, ends);

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

  Pwm matrix("/Users/al/Downloads/pwms/AHR_si.pat", score_threshold);
//  Pwm matrix("/Users/al/Downloads/pwms/TAL1_f2.pat", 13.6);

//  std::stack<std::string> patterns = matrix.getWords(8, 10000000);;
//  while(!matrix.hasMoreWords()) patterns = matrix.getWords(8, 10000000);
  std::vector<std::vector<char> > patterns;
//  patterns.push_back("AGATAATCTTTATTGTCC");
//  patterns = matrix.getWords(8, 10000000);

  AhoCorasickPlus atm;
  
//  std::cout << "#patterns fit: " << patterns.size() << std::endl;
  size_t totalWords = 0;
  while (!matrix.hasMoreWords()) {
    AhoCorasickPlus atm;
    {
      std::vector<std::vector<char> > patterns;
      patterns = matrix.getWords(score_threshold, 50000);
      for (size_t i = 0; i < patterns.size(); i++)
      {
          AhoCorasickPlus::EnumReturnStatus status;
          AhoCorasickPlus::PatternId patId = totalWords++;
          status = atm.addPattern(patterns[i], patId);
          if (status!=AhoCorasickPlus::RETURNSTATUS_SUCCESS)
          {
              std::cerr << "Failed to add: " << std::endl;
          }
      }

      patterns.clear();
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
      
      std::stringstream temp_ss_stream;
      temp_ss_stream << result_folder << sequences.getResultFilename(s_i) << "-temp-" << s_i << time(NULL);
      std::string tempfilename;
      temp_ss_stream >> tempfilename;
      std::ofstream fout(tempfilename);
      
      atm.search(seq, length, false);
      while (atm.findNext(aMatch))
      {
    //    std::cout << "@" << aMatch.position << "\t#" << aMatch.id << "\t" << sample_patterns[aMatch.id] << std::endl;
        fout << aMatch.position << std::endl;
        occurances++;
      }
      
    #ifdef DTIME_PROFILING
      tstop = (double)clock()/CLOCKS_PER_SEC;  
      std::cout << "Search performed in " << tstop - tstart << " seconds" << std::endl;
    #endif
      std::cout << "Occurances found: " << occurances << std::endl;
    }
  }
  return 0;
}
