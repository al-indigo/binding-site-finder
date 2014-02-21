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
#include "utils/prepare_filename.h"



int main(int argc, char **argv) {
//NOTE: results_path should exist (I'm ruling all the filenames stuff with python wrapper - it's easier)
  std::string result_folder("/Users/al/Programming/binding-site-finder/build/results/");
  
  /* 168bytes on model with length 7 */
  /* 205bytes on model with length 9 */
  /* 585bytes on model with length 14 */
  /* 640bytes on model with length 19*/
  /* 870bytes on model with length 25*/ 
  // worst measured case is 42 bytes per character in word.
  

  double score_threshold = -8.0;

  Pwm matrix("/Users/al/Downloads/pwms/AHR_si.pat", score_threshold);
//  Pwm matrix("/Users/al/Downloads/pwms/SOX17_f2.pat", score_threshold);
//  Pwm matrix("/Users/al/Downloads/pwms/TAL1_f2.pat", score_threshold);
//  Pwm matrix("/Users/al/Downloads/pwms/EOMES_f1.pat", score_threshold);
//  Pwm matrix("/Users/al/Downloads/pwms/ZEB1_do.pat", score_threshold);
    
//  Pwm matrix("/Users/al/Downloads/pwms/PAX2_f1.pat", score_threshold);
  
  
  size_t mem_allowed = 8*512; //memory to spend in MB. Don't set it less than READ_BLOCK_SIZE/(1024*1024)
  
  // This is counted in experimental way: it's clear that consumption depends on length, but the coefficient
  // is just experimental for 5 matrices. Real consumption is always less than here (aho-corasick is quite tricky).
  //
  unsigned int patterns_allowed = ( (mem_allowed)*1024*1024 - READ_BLOCK_SIZE) / (32*(matrix.getLength()));
 
  size_t test_start = 0;
  size_t test_end = 249210491;
  
  
  
  std::vector<std::string> filenames;
  std::vector<std::string> descriptions;
  std::vector<std::string> result_filenames;
  std::vector<size_t> starts;
  std::vector<size_t> ends;
  
  filenames.push_back(std::string("/Users/al/Downloads/genome/chr1_upper.plain"));
  descriptions.push_back(std::string("chr1"));
  starts.push_back(test_start);
  ends.push_back(test_end);
  
//  filenames.push_back(std::string("/Users/al/Downloads/genome/chr1_upper.plain"));
//  descriptions.push_back(std::string("chr1 dup"));
//  starts.push_back(test_start);
//  ends.push_back(test_end);
  
  result_filenames.push_back(std::string("chr1-result"));
//  result_filenames.push_back(std::string("chr2-result"));
  
  ChromoVector sequences(filenames, descriptions, result_filenames, starts, ends);
  
  std::vector<std::vector<std::string> > files_to_merge(sequences.size());
//  files_to_merge.reserve(sequences.size());

  /*
  try {} catch (std::exception e) {
#ifdef DDEBUG_PRINT
    std::cerr << e.what();
#endif
  }
*/
  
  std::vector<std::vector<char> > patterns;

  AhoCorasickPlus atm;
  
  size_t totalWords = 0;
  while (!matrix.hasMoreWords()) {
    AhoCorasickPlus atm;
    {
      std::vector<std::vector<char> > patterns;
      patterns = matrix.getWords(score_threshold, patterns_allowed);
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
      for (size_t p_i = 0; p_i < sequences.getNumberOfParts(s_i); p_i++) {
          char * fbuf = new char[FSTREAM_BUF_SIZE];
        
          char * seq = sequences.getSeq(s_i, p_i);
          size_t length = sequences.getPartLength(s_i, p_i);
          size_t offset = sequences.getAbsoluteOffset(s_i, p_i);

          AhoCorasickPlus::Match aMatch;
          size_t occurances = 0;
          std::cout << "Results for " << sequences.getFilename(s_i) << " as " << sequences.getDescription(s_i) << " part #" << p_i << std::endl;

          std::string tempfilename = prepare_filename(result_folder + sequences.getResultFilename(s_i), "-part-", (s_i+1)*10000000 + p_i);
          std::ofstream fout(tempfilename.c_str());
          fout.rdbuf()->pubsetbuf(fbuf, FSTREAM_BUF_SIZE);
          
          files_to_merge[s_i].push_back(tempfilename);

          atm.search(seq, length, false);
          while (atm.findNext(aMatch))
          {
        //    std::cout << "@" << aMatch.position << "\t#" << aMatch.id << "\t" << sample_patterns[aMatch.id] << std::endl;
            fout << aMatch.position + offset - matrix.getLength() << std::endl;
            occurances++;
          }
          
          std::cout << "Occurances found: " << occurances << std::endl;
          fout.close();
          delete [] fbuf;
          sequences.releasePart(s_i, p_i);
      }
    }
  }
  
  std::cout << "Preparing to merge" << std::endl;
  
//   for (int i = 0; i < files_to_merge.size(); i++) {
//     time_t t = time(NULL);
//     std::string tempfilename = prepare_filename(result_folder + result_filenames[i], std::string("-temp-merge-"), 0, &t);
//     std::ofstream fout(tempfilename.c_str());
//     fout.close();
//     
//     for (int j = 0; j < files_to_merge[i].size(); j++) {
//       std::string read_file = prepare_filename(result_folder + result_filenames[i], std::string("-temp-merge-"), j, &t);
//       std::string write_file = prepare_filename(result_folder + result_filenames[i], std::string("-temp-merge-"), j+1, &t);
//       
//       merge_sort(files_to_merge[i][j], read_file, write_file);
//     }
//     
//   }
  
  for (int i = 0; i < files_to_merge.size(); i++) {
    time_t t = time(NULL);
    std::string tempfilename = prepare_filename(result_folder + result_filenames[i], std::string("-temp-merge-stl-"), 0, &t);
    std::ofstream fout(tempfilename.c_str());
    fout.close();
    
    for (int j = 0; j < files_to_merge[i].size(); j++) {
      std::string read_file = prepare_filename(result_folder + result_filenames[i], std::string("-temp-merge-stl-"), j, &t);
      std::string write_file = prepare_filename(result_folder + result_filenames[i], std::string("-temp-merge-stl-"), j+1, &t);
      
      merge_sort(files_to_merge[i][j], read_file, write_file);
    }
    
  }
  
  {
    // read part of merged file ints
    // read corresponding words
    // compute p-values and scores for these words
    // write this bunch of values to result file
  }
  
  
  return 0;
}
