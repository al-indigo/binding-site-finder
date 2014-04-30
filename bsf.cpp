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
#include <sstream>
#include <thread>
#include <functional>
#include <mutex>

#include "ahoc/AhoCorasickPlus.h"
#include "structs/pwm.h"
#include "structs/chromovector.h"
#include "utils/mergesort_infile.h"
#include "utils/prepare_filename.h"
#include "utils/jsonxx.h"
#include "functions/functions.h"

enum methods {naive, ahoc};

/* IMPORTANT NOTE: 214600 + 0.0006868 is an experimental value for out I/O system in production. In your situation it may differ (a lot!) especially if you use SSD.
 * How to measure it:
 *  1) choose a model (with size 13 is the optimal one because nothing else iterefers the results (smaller give more results, larger take too much time for matrix traversal)
 *  2) measure time for running this model with p-values 0.00001 up to 0.001 with both methods, remember also number of words that you are looking for (it's p-value*voc_volume)
 *  3) make a plot for dependency time/number of words. 
 *  4) fit a trace for convinience -- it linear for both cases
 *  5) intersection point is the value where algorithms work the same time
 *  6) make steps 2-5 for another chromosome with different size (if you choose chromosomes with size difference > 2, it would be better. Also try to choose chromos with size > 96MB (it's buffer by default))
 *  7) fit a trace for this two points by (number of words to find/intersection point).
 *  8) coefficients for this linear are your calibration value. */
#define CALIBRATION_VALUE(total_length) 214600 + 0.0006868 * total_length

void naive_walk(size_t& s_i, size_t& start, size_t&stop, std::vector<double>& scores, std::vector<bool>& strand, std::vector<double>& pvalues, std::vector<size_t>& matched, Pwm& matrix, ChromoVector& sequences, size_t& offset) {

  sequences.getWordScoresVector(s_i, offset + start, offset + stop, matrix, scores, strand, matched);

  pvalues.resize(scores.size());
  
  matrix.getPValuesPlain(scores, pvalues);
}

int predict(size_t mem_allowed,
            std::string matrix_filename, 
            const double p_value, 
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
  
  matrix_optimization_type type;
  if (('A') != 0x41 ||
      ('C') != 0x43 ||
      ('G') != 0x47 ||
      ('T') != 0x54 ||
      ('N') != 0x4E ||
      ('a'^0x20) != 0x41 ||
      ('c'^0x20) != 0x43 ||
      ('g'^0x20) != 0x47 ||
      ('t'^0x20) != 0x54 ||
      ('n'^0x20) != 0x4E) {
    std::cout << "WARNING! Your ASCII symbol codes are not standard; you must be sure that your genome files are uppercased otherwise results may be incorrect\n";
    type = pat; 
  } else {
    type = pat_bit_optimization;
  }

  Pwm matrix(matrix_filename.c_str(), p_value, type);
  
  size_t total_length = 0;
  for (auto& i: starts) {
    total_length += ends[i] - starts[i];
  }
  
  methods method;

  if (matrix.getLength() > 19 || matrix.getNumberOfWords() > CALIBRATION_VALUE(total_length)) {
    method = naive;
  } else {
    method = ahoc;
  }
 
  // This is counted in experimental way: it's clear that consumption depends on length, but the coefficient
  // is just experimental for 5 matrices. Real consumption is always less than here (aho-corasick is quite tricky).
  unsigned int patterns_allowed = ( mem_allowed - READ_BLOCK_SIZE) / (32*(matrix.getLength()));
  
  ChromoVector sequences(filenames, chromonames, starts, ends, matrix.getLength());

  std::vector<std::vector<std::string> > files_to_merge(sequences.size());
  method = naive;
  if (method == naive) {
    for (size_t s_i = 0; s_i < sequences.size(); s_i++) {
      for (size_t p_i = 0; p_i < sequences.getNumberOfParts(s_i); p_i++) {
        sequences.getSeq(s_i, p_i);
        size_t offset = sequences.getAbsoluteOffset(s_i, p_i);
        int numthreads = 4;
        std::vector<std::vector<double> > scores(numthreads), pvalues(numthreads);
        std::vector<std::vector<size_t> > matched(numthreads);
        std::vector<std::vector<bool> > strand(numthreads);
        std::vector<size_t> st(numthreads), en(numthreads);
        for (int i = 0; i < numthreads; i++) {
          st[i] = i*(sequences.getPartLength(s_i, p_i)/numthreads);
          en[i] = (i+1)*(sequences.getPartLength(s_i, p_i)/numthreads);
        }
        en[numthreads-1] = sequences.getPartLength(s_i, p_i) - matrix.getLength();
        
        std::vector<std::thread> threads;
        for (int i = 0; i < numthreads; i++) {
          threads.push_back(std::thread(naive_walk, std::ref(s_i), std::ref(st[i]), std::ref(en[i]), std::ref(scores[i]), std::ref(strand[i]), std::ref(pvalues[i]), std::ref(matched[i]), std::ref(matrix), std::ref(sequences), std::ref(offset)));
        }
      
        FILE * fout = fopen((result_folder + result_filename).c_str(), "a");
        for (int i = 0; i < numthreads; i++) {
          threads[i].join();
          format_bed(fout, chromonames[s_i], matched[i], pvalues[i], scores[i], strand[i], st[i]);
          write_status(100.0 * (((double)(s_i + 1) * (p_i + 1)) / ((double) sequences.size() * sequences.getNumberOfParts(s_i))), status_folder, status_filename, "searching (naive)", "" );
        }
        fflush(fout);
        fclose(fout);
      }
    }
    write_status(100.0, status_folder, status_filename, "task complete", (std::string("http://bsf.at.ispras.ru/result-files/") + result_filename).c_str());
  } else {
    size_t total_words = 0;
    
    while (matrix.hasMoreWords()) {
      AhoCorasickPlus atm;
      {
        init_ahoc(atm, matrix, patterns_allowed, total_words);
      }
      
      atm.finalize();
      //TODO: try to make correct copy constructor for finalized automata (to make parallel searches)
      for (size_t s_i = 0; s_i < sequences.size(); s_i++) {
        for (size_t p_i = 0; p_i < sequences.getNumberOfParts(s_i); p_i++) {
          ahoc_search(s_i, p_i, atm, matrix, sequences, result_folder, result_filename, total_words, files_to_merge);
          write_status(60.0 * ((double)(matrix.getNumberOfWordsFound() * (s_i + 1) * (p_i + 1)) / (double)(matrix.getNumberOfWords() * sequences.size() * sequences.getNumberOfParts(s_i))), status_folder, status_filename, "searching", "");
        }
      }
    }
  
    std::vector<std::string> merged_files;
    
    for (size_t i = 0; i < files_to_merge.size(); i++) {
      merge_files(i, result_folder, result_filename, files_to_merge, merged_files);
      
      write_status(60.0 + 15.0 * (double)i/(double)files_to_merge.size(), status_folder, status_filename, "merging files", "");
    }
 
    FILE * fout = fopen((result_folder + result_filename).c_str(), "a");
    for (size_t i = 0; i < merged_files.size(); i++) {
      FILE * fin = fopen(merged_files[i].c_str(), "r");
      while (!feof(fin)) {
        std::vector<bool> strand;
        std::vector<size_t> matched;
        std::vector<double> scoresFw, scoresRev, pvaluesFw, pvaluesRev;
            
//        recalc_scores_p_values(fin, positions_to_read_again, mem_allowed, matrix, sequences, i, pvaluesFw, pvaluesRev, scoresFw, scoresRev);
        recalc_scores_p_values(fin, mem_allowed, matrix, sequences, i, pvaluesFw, scoresFw, strand, matched);
        if (matched.size() == 0) continue;
              
        format_bed(fout, chromonames[i], matched, pvaluesFw, scoresFw, strand, size_t(0));
        
      }
      
      fclose(fin);
      remove(merged_files[i].c_str());
      
      write_status(75.0 + 25.0 * (double) i / (double) merged_files.size(), status_folder, status_filename, "recomputing scores and p-values", "");
    }
    fflush(fout);
    fclose(fout);

    write_status(100.0, status_folder, status_filename, "task complete", (std::string("http://bsf.at.ispras.ru/result-files/") + result_filename).c_str());
  }
  return 0;
}


double getThres(std::string matrix_filename, double p_value) {
  Pwm matrix(matrix_filename, p_value);
  std::cout << "Threshold " << matrix.getThreshold() << " words approximately " << matrix.getNumberOfWordLeft() << std::endl;

  return matrix.getThreshold();
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
  
  getThres(matrix_filename, p_value);
  
  predict(mem_allowed, matrix_filename, p_value, result_folder, result_filename, status_folder, status_filename, filenames, chromonames, starts, ends);

  return 0;
}
