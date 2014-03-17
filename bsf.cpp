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
#include <cstdio>

#include "ahoc/AhoCorasickPlus.h"
#include "structs/pwm.h"
#include "structs/chromovector.h"
#include "utils/mergesort_infile.h"
#include "utils/prepare_filename.h"
#include "utils/jsonxx.h"

int predict(size_t mem_allowed,
            std::string matrix_filename, 
            double p_value, 
            std::string result_folder,
            std::string result_filename,
            std::string status_folder,
            std::string status_filename,
            std::vector<std::string> filenames, 
            std::vector<std::string> chromonames, 
            std::vector<size_t> starts, 
            std::vector<size_t> ends) {
  
//NOTE: results_path should exist (I'm ruling all the filenames stuff with python wrapper - it's easier)
 
  Pwm matrix(matrix_filename.c_str(), p_value);
 
  // This is counted in experimental way: it's clear that consumption depends on length, but the coefficient
  // is just experimental for 5 matrices. Real consumption is always less than here (aho-corasick is quite tricky).
  //
  unsigned int patterns_allowed = ( mem_allowed - READ_BLOCK_SIZE) / (32*(matrix.getLength()));
  
  ChromoVector sequences(filenames, chromonames, starts, ends);
  
  std::vector<std::vector<std::string> > files_to_merge(sequences.size());

  std::vector<std::vector<char> > patterns;

  AhoCorasickPlus atm;
  
  size_t totalWords = 0;
  double voc_volume = (double)(matrix.getLength()^4) * p_value;
  double status = 0.0;
  
  while (!matrix.hasMoreWords()) {
    AhoCorasickPlus atm;
    {
      std::vector<std::vector<char> > patterns;
      patterns = matrix.getWords(patterns_allowed);
      for (size_t i = 0; i < patterns.size(); i++)
      {
          AhoCorasickPlus::EnumReturnStatus res;
          AhoCorasickPlus::PatternId patId = totalWords++;
          res = atm.addPattern(patterns[i], patId);
          if (res!=AhoCorasickPlus::RETURNSTATUS_SUCCESS)
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
            result_vector.push_back(aMatch.position + offset - matrix.getLength());
            occurances++;
          }
          
          std::cout << "Occurances found: " << occurances << std::endl;
          if (occurances > 0) {
            std::string tempfilename = prepare_filename(result_folder + result_filename, "bsf-part-", (s_i+1)*10000000 + p_i);
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
          status = 60 * (double)totalWords / voc_volume;
          std::ofstream status_out((status_folder + status_filename).c_str());
          status_out << "{\"percent done\": " << (int)status << ", \"result\": [\"\"], \"explain\": \"searching\"}";
          status_out.close();
      }
    }
  }
  
  std::vector<std::string> merged_files;
  
  for (int i = 0; i < files_to_merge.size(); i++) {
    time_t t = time(NULL);
    std::string tempfilename = prepare_filename(result_folder + result_filename, std::string("-temp-merge-stl-"), 10000000*(i+1) + 0, &t);
    std::ofstream fout(tempfilename.c_str());
    fout.close();
    
    for (int j = 0; j < files_to_merge[i].size(); j++) {
      std::string read_file = prepare_filename(result_folder + result_filename, std::string("-temp-merge-stl-"), 10000000*(i+1) + j, &t);
      std::string write_file = prepare_filename(result_folder + result_filename, std::string("-temp-merge-stl-"), 10000000*(i+1) + j+1, &t);
      
      merge_sort(files_to_merge[i][j], read_file, write_file);
    }
    merged_files.push_back(prepare_filename(result_folder + result_filename, std::string("-temp-merge-stl-"), 10000000*(i+1) + files_to_merge[i].size(), &t));
    
    status = 60.0 + 15.0 * (double)i/(double)files_to_merge.size();
    std::ofstream status_out((status_folder + status_filename).c_str());
    status_out << "{\"percent done\": " << (int)status << ", \"result\": [\"\"], \"explain\": \"merging files\"}";
    status_out.close();
  }

  std::ofstream fout((result_folder + result_filename).c_str(), std::fstream::app);
  for (int i = 0; i < merged_files.size(); i++) {
    std::ifstream fin(merged_files[i].c_str());
    fin >> std::setbase(16);
    while (fin) {
      std::set<size_t> positions_to_read_again;
      size_t buf = 0;
      for (int k = 0; k < mem_allowed * 1024 * 1024 / 8 && fin >> buf; k++) {
        positions_to_read_again.insert(buf);
      }
      if (positions_to_read_again.empty()) {
        continue;
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
      std::cout << "Result file: " << result_filename << std::endl;
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
        
        fout      << chromonames[i] << "\t" 
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
    remove(merged_files[i].c_str());
    
    status = 75.0 + 25.0 * (double) i / (double) merged_files.size();
    std::ofstream status_out((status_folder + status_filename).c_str());
    status_out << "{\"percent done\": " << (int)status << ", \"result\": [\"\"], \"explain\": \"recomputing scores and p-values\"}";
    status_out.close();    
  }
  fout.close();
  
  std::ofstream status_out((status_folder + status_filename).c_str());
  status_out << "{\"percent done\": 100, \"result\":" << "[\"http://bsf.at.ispras.ru/result-files/" << result_filename << "\"], \"explain\": \"task complete\"}";
  status_out.close();    
 
  return 0;
}



int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "You must provide a json file with task as an input. You may look at the example.json distributed with sources" << std::endl;
    exit(1);
  }
  
  std::string cli(argv[1]);
  
  std::ifstream options(cli.c_str());
  
  jsonxx::Object task;
  assert(task.parse(options));
  
  size_t mem_allowed = task.get<jsonxx::Number>("memory");
  std::string matrix_filename = task.get<jsonxx::String>("matrix");
  double p_value = task.get<jsonxx::Number>("p-value");
  std::string result_folder = task.get<jsonxx::String>("result-folder");
  std::string result_filename = task.get<jsonxx::String>("result-filename");
  std::string status_folder = task.get<jsonxx::String>("status-folder");
  std::string status_filename = task.get<jsonxx::String>("status-filename");
  
  std::vector<std::string> filenames;
  std::vector<std::string> chromonames;
  std::vector<size_t> starts;
  std::vector<size_t> ends;
  
  if (task.get<jsonxx::Array>("tasks").size() == 0) {
    std::cerr << "Provided Json file doesn't contain tasks" << std::endl;
    exit(1);
  }
  
  for (size_t i = 0; i < task.get<jsonxx::Array>("tasks").size(); i++) {
    filenames.push_back(task.get<jsonxx::Array>("tasks").get<jsonxx::Object>(i).get<jsonxx::String>("chromofile"));
    chromonames.push_back(task.get<jsonxx::Array>("tasks").get<jsonxx::Object>(i).get<jsonxx::String>("chromoname"));
    starts.push_back(task.get<jsonxx::Array>("tasks").get<jsonxx::Object>(i).get<jsonxx::Number>("start"));
    ends.push_back(task.get<jsonxx::Array>("tasks").get<jsonxx::Object>(i).get<jsonxx::Number>("end"));
  }
  
  return predict(mem_allowed, matrix_filename, p_value, result_folder, result_filename, status_folder, status_filename, filenames, chromonames, starts, ends);
}