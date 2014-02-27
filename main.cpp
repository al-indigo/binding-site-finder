#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

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
  

//  double score_threshold = -8.0;
//  double p_value = 0.000000001;
    double p_value = 0.00107763671875;

//  Pwm matrix("/Users/al/Programming/perfectosape/test_data/pwm/KLF4_f2.pwm", p_value);
  
  Pwm matrix("/Users/al/Downloads/pwms/AHR_si.pat", p_value);
//  Pwm matrix("/Users/al/Downloads/pwms/SOX17_f2.pat", score_threshold);
//  Pwm matrix("/Users/al/Downloads/pwms/TAL1_f2.pat", p_value);
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
  
//   filenames.push_back(std::string("/Users/al/Downloads/genome/chr1_upper.plain"));
//   descriptions.push_back(std::string("chr1 dup"));
//   starts.push_back(test_start);
//   ends.push_back(test_end);
  
  result_filenames.push_back(std::string("chr1-result"));
//   result_filenames.push_back(std::string("chr2-result"));
  
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
      patterns = matrix.getWords(patterns_allowed);
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
          
          std::vector<size_t> result_vector;
          result_vector.reserve(READ_BLOCK_SIZE);
        
          char * seq = sequences.getSeq(s_i, p_i);
          size_t length = sequences.getPartLength(s_i, p_i);
          size_t offset = sequences.getAbsoluteOffset(s_i, p_i);

          AhoCorasickPlus::Match aMatch;
          size_t occurances = 0;
          std::cout << "Results for " << sequences.getFilename(s_i) << " as " << sequences.getDescription(s_i) << " part #" << p_i << std::endl;

          atm.search(seq, length, false);
          while (atm.findNext(aMatch))
          {
        //    std::cout << "@" << aMatch.position << "\t#" << aMatch.id << "\t" << sample_patterns[aMatch.id] << std::endl;
            result_vector.push_back(aMatch.position + offset - matrix.getLength());
            occurances++;
          }
          
          std::cout << "Occurances found: " << occurances << std::endl;
          if (occurances > 0) {
            std::string tempfilename = prepare_filename(result_folder + sequences.getResultFilename(s_i), "-part-", (s_i+1)*10000000 + p_i);
            std::ofstream fout(tempfilename.c_str());
            char * fbuf = new char[FSTREAM_BUF_SIZE];
            fout.rdbuf()->pubsetbuf(fbuf, FSTREAM_BUF_SIZE);
            fout << std::setbase(16);

            files_to_merge[s_i].push_back(tempfilename);
            
            for (std::vector<size_t>::iterator z = result_vector.begin(); z != result_vector.end(); ++z) {
              fout << *z << std::endl;
            }
            
            fout.close();
            delete [] fbuf;
          }
          
          sequences.releasePart(s_i, p_i);
      }
    }
  }
  
  std::cout << "Preparing to merge" << std::endl;

  std::vector<std::string> merged_files;
  
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
    merged_files.push_back(prepare_filename(result_folder + result_filenames[i], std::string("-temp-merge-stl-"), files_to_merge[i].size(), &t));
  }
  
  for (int i = 0; i < merged_files.size(); i++) {
    std::ifstream fin(merged_files[i]);
    std::ofstream fout(result_filenames[i]);
    fin >> std::setbase(16);
    while (fin) {
      std::set<size_t> positions_to_read_again;
      for (int k = 0; k < mem_allowed * 1024 * 1024 / 8 && fin; k++) {
        size_t buf = 0;
        fin >> buf;
        if (buf == 0 && k != 0) { continue; }
          
        positions_to_read_again.insert(buf);
        
        if (buf == 0) {
          std::cout << "I'm null";
        }
      }
      
      std::vector<std::vector<char> > words_as_paths(positions_to_read_again.size());
      std::vector<double> scoresFw;
      std::vector<double> scoresRev;
      std::vector<double> pvaluesFw;
      std::vector<double> pvaluesRev;
      
      sequences.getWordsAsPaths(i, positions_to_read_again, matrix.getLength(), words_as_paths);
      matrix.getScores(words_as_paths, scoresFw, scoresRev);
      
      pvaluesFw = matrix.getPValues(scoresFw);
      pvaluesRev = matrix.getPValues(scoresRev);
      
      std::set<size_t>::iterator positerator = positions_to_read_again.begin()++;
      std::cout.precision(8);
      std::cout << "Result file: " << result_filenames[i] << std::endl;
      for (int k = 0; k < positions_to_read_again.size() ; k++) {
        double pv = std::min(pvaluesFw[k], pvaluesRev[k]);
        int narrowPeakScore = 0;
        if (pv > 0.001) {
          narrowPeakScore = 100; 
        } else if (pv < 0.00001) {
          narrowPeakScore = 1000;
        } else {
          narrowPeakScore = 100 + (int) (900 * ((0.001 - pv)/(0.00099) ) );
        }
        
        std::cout << descriptions[i] << "\t" 
                  << *positerator << "\t" << *positerator << "\t"
                  << ".\t"
                  << narrowPeakScore << "\t" 
                  << (scoresFw[k] > scoresRev[k] ? "+" : "-") << "\t" 
                  << std::max(scoresFw[k], scoresRev[k]) << "\t" 
                  << pv << "\t" << "-1\t-1" << std::endl;
        positerator++;
      }
      
      positions_to_read_again.clear();
    }
    
    fin.close();
  }
  
  {
    // read part of merged file ints
    // read corresponding words
    // compute p-values and scores for these words
    // write this bunch of values to result file
  }
  
  
  return 0;
}
