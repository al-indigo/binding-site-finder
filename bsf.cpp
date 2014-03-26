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
#include "functions/functions.h"

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
  {
    std::ofstream fout((result_folder + result_filename).c_str(), std::fstream::out);
    fout.close();
  }
 
  Pwm matrix(matrix_filename.c_str(), p_value);
 
  // This is counted in experimental way: it's clear that consumption depends on length, but the coefficient
  // is just experimental for 5 matrices. Real consumption is always less than here (aho-corasick is quite tricky).
  //
  unsigned int patterns_allowed = ( mem_allowed - READ_BLOCK_SIZE) / (32*(matrix.getLength()));
  
  ChromoVector sequences(filenames, chromonames, starts, ends);
  
  std::vector<std::vector<std::string> > files_to_merge(sequences.size());

  AhoCorasickPlus atm;
  
  size_t total_words = 0;
  double voc_volume = (double)(matrix.getLength()^4) * p_value;
  double status = 0.0;
  
  while (!matrix.hasMoreWords()) {
    AhoCorasickPlus atm;
    {
      init_ahoc(atm, matrix, patterns_allowed, total_words);
    }
    
    atm.finalize();
    
    for (size_t s_i = 0; s_i < sequences.size(); s_i++) {
      for (size_t p_i = 0; p_i < sequences.getNumberOfParts(s_i); p_i++) {
        ahoc_search(s_i, p_i, atm, matrix, sequences, result_folder, result_filename, total_words, files_to_merge);
        write_status(60 * (double)total_words / voc_volume, status_folder, status_filename, "searching", "");
      }
    }
  }
  
  std::vector<std::string> merged_files;
  
  for (int i = 0; i < files_to_merge.size(); i++) {
    merge_files(i, result_folder, result_filename, files_to_merge, merged_files);
    
    write_status(60.0 + 15.0 * (double)i/(double)files_to_merge.size(), status_folder, status_filename, "merging files", "");
  }

  std::ofstream fout((result_folder + result_filename).c_str(), std::fstream::app);
  for (int i = 0; i < merged_files.size(); i++) {
    std::ifstream fin(merged_files[i].c_str());
    fin >> std::setbase(16);
    while (fin) {
      std::set<size_t> positions_to_read_again;
      std::vector<double> scoresFw, scoresRev, pvaluesFw, pvaluesRev;
          
      recalc_scores_p_values(fin, positions_to_read_again, mem_allowed, matrix, sequences, i, pvaluesFw, pvaluesRev, scoresFw, scoresRev);
      if (positions_to_read_again.size() == 0) continue;
      
      format_bed(fout, chromonames[i], positions_to_read_again, pvaluesFw, pvaluesRev, scoresFw, scoresRev);
    }
    
    fin.close();
    remove(merged_files[i].c_str());
    
    write_status(75.0 + 25.0 * (double) i / (double) merged_files.size(), status_folder, status_filename, "recomputing scores and p-values", "");
  }
  fout.close();
  
  write_status(100.0, status_folder, status_filename, "task complete", (std::string("http://bsf.at.ispras.ru/result-files/") + result_filename).c_str());
 
  return 0;
}



int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "You must provide a json file with task as an input. You may look at the example.json distributed with sources" << std::endl;
    exit(1);
  }
  
  std::string cli(argv[1]);
  
  std::ifstream options(cli.c_str());
  
  if (!options) {
    exit(1);
  }
  
  jsonxx::Object task;
  task.parse(options);
  
  size_t mem_allowed = 1024*1024*task.get<jsonxx::Number>("memory");
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