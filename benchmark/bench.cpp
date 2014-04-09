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

#include "../ahoc/AhoCorasickPlus.h"
#include "../structs/pwm.h"
#include "../structs/chromovector.h"
#include "../utils/mergesort_infile.h"
#include "../utils/prepare_filename.h"
#include "../functions/functions.h"
#include "../utils/jsonxx.h"

#include "hayai/hayai.hpp"

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif


void parse_options(const char * task_filename, 
                   size_t& mem_allowed, 
                   size_t& genome_len,
                   size_t& scores_to_test,
                   size_t& p_values_to_test,
                   std::string& matrix_filename, 
                   double& p_value) {
  
  std::ifstream options(task_filename);
  if (!options) {
    std::cerr << "Task doesn't exist: " << task_filename << std::endl;
    exit(1);
  }
  
  jsonxx::Object task;
  task.parse(options);
  
  mem_allowed = 1024*1024*task.get<jsonxx::Number>("memory");
  genome_len = task.get<jsonxx::Number>("genome-len");
  matrix_filename = task.get<jsonxx::String>("matrix");
  p_value = task.get<jsonxx::Number>("p-value");
  scores_to_test = task.get<jsonxx::Number>("scores-to-test");
  p_values_to_test = task.get<jsonxx::Number>("p-values-to-test");
}


std::string get_working_path()
{
   char path[1024];   
   char real_path[1024];
#ifdef __APPLE__
   uint32_t size = sizeof(path);
   if (_NSGetExecutablePath(path, &size) == 0) {
     realpath(path, real_path);
     return std::string(real_path);
   } else {
     exit(-1);
   }
#endif
   
#ifdef __linux
   return std::string("linux");
#endif
   
   return std::string("/dev/null/");
}

size_t check_file_exists(std::string filename) {
  std::ifstream f(filename.c_str());
  if (!f.good()) {
    f.close();
    return 0;
  } else {
    size_t counter = 0;
    char temp;
    while(f >> temp) {
      counter++;
    }
    f.close();
    return counter;
  }
}


struct GlobalParameters {
  size_t mem_allowed;
  unsigned int patterns_allowed;
  size_t genome_len;
  std::string matrix_filename;
  double p_value;
  size_t scores_to_test;
  size_t p_values_to_test;
  std::string result_folder;
  std::string result_filename;
  std::string status_folder;
  std::string status_filename;
  
  std::vector<std::string> filenames;
  std::vector<std::string> chromonames;
  std::vector<size_t> starts;
  std::vector<size_t> ends;
  
  GlobalParameters() {
      parse_options("benchmark.json", mem_allowed, genome_len, scores_to_test, p_values_to_test, matrix_filename, p_value);
      
      std::cout << hayai::Console::TextGreen << "[Testing parameters: ]" << hayai::Console::TextDefault << std::endl;
      std::cout << "\tMemory allowed to use: " << mem_allowed << " bytes" << std::endl;
      
      Pwm matrix(matrix_filename, p_value);
      patterns_allowed = ( mem_allowed - READ_BLOCK_SIZE) / (32*(matrix.getLength()));
      
      std::cout << "\tThis means that we will process maximum " << patterns_allowed << " patterns in each point of time." << std::endl;
      std::cout << "\tModel is: " << matrix_filename << std::endl;
      std::cout << "\tSelected p_value: " << p_value << std::endl;
      
      size_t len = check_file_exists(get_working_path() + std::string("-chr_random.plain"));
      if (len == 0) {
        std::cout << "\tTest genome missing; generating random genome with length: " << genome_len << std::endl;
        std::string chromofile = get_working_path() + std::string("-chr_random.plain");
        std::ofstream chromo(chromofile.c_str());
        static const char alphabet[] = "ACGT";
        for (size_t i=0; i < genome_len; i++) {
          chromo << alphabet[rand() % (sizeof(alphabet)-1)];
        }
        chromo.close();
      } else {
        std::cout << "\tTest genome is present; it's length equals " << len << std::endl;
        genome_len = len;
      }
      
      if (scores_to_test > genome_len) scores_to_test = genome_len - 100;
      
      std::cout << hayai::Console::TextGreen << "[Parameters for benchmark initialized successfully]" << hayai::Console::TextDefault << std::endl;
      
  }
  
};

struct GlobalParameters global_params;



class GetWordsTest : public hayai::Fixture {
public:
  virtual void SetUp() {
    this->matrix = new Pwm(global_params.matrix_filename, global_params.p_value);
  }
  virtual void TearDown()
  {
    delete this->matrix;
  }
  
  void getWords(optimization_type type) {
    if (type != classic && matrix->getLength() > 18) {
      std::cout << "This model is too big, the process will never end, stopping this test" << std::endl;
      return;
    }

    while (matrix->hasMoreWords()) {
      (matrix->getWords(global_params.patterns_allowed, type)).size();
    }
    
  }
  
  Pwm* matrix;
};


class InitAhoCAutomata : public hayai::Fixture {
public:
  virtual void SetUp() {
    size_t total_words = 0;
    
    this->matrix = new Pwm(global_params.matrix_filename, global_params.p_value);

    this->atm = new AhoCorasickPlus();
    
    init_ahoc(*(this->atm), *(this->matrix), global_params.patterns_allowed, total_words);
  }
  
  virtual void TearDown()
  {
    delete this->matrix;
    delete this->atm;
  }
  
  void initAutomata() {
    this->atm->finalize();
  }
  
  Pwm* matrix;
  AhoCorasickPlus* atm;
};

class SearchAhoC : public hayai::Fixture {
public:
  virtual void SetUp() {
    this->matrix = new Pwm(global_params.matrix_filename, global_params.p_value);
    this->atm = new AhoCorasickPlus();
    init_ahoc(*(this->atm), *(this->matrix), global_params.patterns_allowed, total_words);
    atm->finalize();
    
    std::vector<string> filenames, chromonames;
    std::vector<size_t> starts, ends;
    filenames.push_back(get_working_path() + std::string("-chr_random.plain"));
    chromonames.push_back(std::string("bench_chromo"));
    starts.push_back(0);
    ends.push_back(global_params.genome_len-1);
    
    sequences = new ChromoVector(filenames, chromonames, starts, ends);
    total_words=0;
    files_to_merge = new std::vector<std::vector<std::string> >;
  }
  
  virtual void TearDown() {
    
    for (int i = 0; i < this->files_to_merge->size(); i++) {
      for (int j = 0; j < this->files_to_merge[i].size(); j++) {
        remove((*(this->files_to_merge))[i][j].c_str());
      }
    }
    
    delete files_to_merge;
    delete sequences;
    delete matrix;
    delete atm;
  }
  Pwm* matrix;
  AhoCorasickPlus* atm;
  size_t total_words;
  ChromoVector* sequences;
  std::vector<std::vector<std::string> >* files_to_merge;
};




class RecalcScores : public hayai::Fixture {
public:
  virtual void SetUp() {
    this->matrix = new Pwm(global_params.matrix_filename, global_params.p_value);
    
    std::vector<string> filenames, chromonames;
    std::vector<size_t> starts, ends;
    filenames.push_back(get_working_path() + std::string("-chr_random.plain"));
    chromonames.push_back(std::string("bench_chromo"));
    starts.push_back(0);
    ends.push_back(global_params.genome_len-1);
    
    sequences = new ChromoVector(filenames, chromonames, starts, ends);
    for (size_t i = 0; i < global_params.scores_to_test - this->matrix->getLength() -1 ; i++) {
      stream << i << std::endl;
    }
    
    size_t buf = 0;
    for (int k = 0; k < global_params.mem_allowed / 8 && stream >> buf; k++) {
      this->positions_to_read_again.push_back(buf);
    }
    if (positions_to_read_again.empty()) {
      std::cout << "ERROR!" << std::endl;
    }
    
    words_as_paths.resize(global_params.scores_to_test - this->matrix->getLength());
    this->sequences->getWordsAsPaths(0, positions_to_read_again, this->matrix->getLength(), words_as_paths);
  }
  
  virtual void TearDown() {
    words_as_paths.clear();
    this->positions_to_read_again.clear();
    this->pvaluesFw.clear();
    this->pvaluesRev.clear();
    this->scoresFw.clear();
    this->scoresRev.clear();
    delete sequences;
    delete matrix;
  }
  
  void recalc_score(optimization_type type) {
    std::cout << hayai::Console::TextCyan << "Calcucating " << positions_to_read_again.size() << " scores, word length is " << matrix->getLength() << hayai::Console::TextDefault << std::endl;
    this->matrix->getScores(this->words_as_paths, scoresFw, scoresRev, type);
  }
  
  void recalc_full() {
     recalc_scores_p_values(this->stream, 
                            this->positions_to_read_again, 
                            global_params.mem_allowed, 
                            *(this->matrix), 
                            *(this->sequences), 
                            0, 
                            this->pvaluesFw, 
                            this->pvaluesRev, 
                            this->scoresFw, 
                            this->scoresRev);
  }
  
  std::vector<double> pvaluesFw, pvaluesRev, scoresFw, scoresRev;
  std::stringstream stream;
  std::vector<size_t> positions_to_read_again;
  Pwm* matrix;
  ChromoVector* sequences;
  std::vector<std::vector<char> > words_as_paths;
};

class ReadWords : public hayai::Fixture {
public:
  virtual void SetUp() {
    this->matrix = new Pwm(global_params.matrix_filename, global_params.p_value);
    
    std::vector<string> filenames, chromonames;
    std::vector<size_t> starts, ends;
    filenames.push_back(get_working_path() + std::string("-chr_random.plain"));
    chromonames.push_back(std::string("bench_chromo"));
    starts.push_back(0);
    ends.push_back(global_params.genome_len-1);
    
    sequences = new ChromoVector(filenames, chromonames, starts, ends);
    for (size_t i = 0; i < global_params.scores_to_test - this->matrix->getLength() -1 ; i++) {
      stream << i << std::endl;
    }

  }
  
  virtual void TearDown() {
    words_as_paths.clear();
    this->positions_to_read_again.clear();
    this->pvaluesFw.clear();
    this->pvaluesRev.clear();
    this->scoresFw.clear();
    this->scoresRev.clear();
    delete sequences;
    delete matrix;
  }
  
  void read_new() {
    size_t buf = 0;
    words_as_paths.resize(global_params.scores_to_test - this->matrix->getLength(), std::vector<char>(this->matrix->getLength()) );
    for (int k = 0; k < global_params.mem_allowed / 8 && stream >> buf; k++) {
      this->sequences->getWordAsPathTest(0, buf, this->matrix->getLength(), temp);
      this->words_as_paths.push_back(temp);
    }
    
    this->matrix->getScores(this->words_as_paths, scoresFw, scoresRev, classic);
  }
  
  void read_old() {
    size_t buf = 0;
    for (int k = 0; k < global_params.mem_allowed / 8 && stream >> buf; k++) {
      this->positions_to_read_again.push_back(buf);
    }
    if (positions_to_read_again.empty()) {
      std::cout << "ERROR!" << std::endl;
    }
    
    words_as_paths.resize(global_params.scores_to_test - this->matrix->getLength(), std::vector<char>(this->matrix->getLength()) );
    this->sequences->getWordsAsPaths(0, positions_to_read_again, this->matrix->getLength(), words_as_paths);
    this->matrix->getScores(this->words_as_paths, scoresFw, scoresRev, classic);
  }
  
  std::vector<char> temp;
  std::vector<double> pvaluesFw, pvaluesRev, scoresFw, scoresRev;
  std::stringstream stream;
  std::vector<size_t> positions_to_read_again;
  Pwm* matrix;
  ChromoVector* sequences;
  std::vector<std::vector<char> > words_as_paths;
};


class CalcPvalues : public hayai::Fixture {
public:
  virtual void SetUp() {
    this->matrix = new Pwm(global_params.matrix_filename, global_params.p_value);
    
    std::vector<string> filenames, chromonames;
    std::vector<size_t> starts, ends;
    filenames.push_back(get_working_path() + std::string("-chr_random.plain"));
    chromonames.push_back(std::string("bench_chromo"));
    starts.push_back(0);
    ends.push_back(global_params.genome_len-1);
    
    sequences = new ChromoVector(filenames, chromonames, starts, ends);
    for (size_t i = 0; i < global_params.p_values_to_test ; i++) {
      stream << i << std::endl;
    }
    
    size_t buf = 0;
    words_as_paths.resize(global_params.p_values_to_test, std::vector<char>(this->matrix->getLength()) );
    for (int k = 0; k < global_params.mem_allowed / 8 && stream >> buf; k++) {
      this->sequences->getWordAsPathTest(0, buf, this->matrix->getLength(), temp);
      this->words_as_paths.push_back(temp);
    }
    
    this->matrix->getScores(this->words_as_paths, scoresFw, scoresRev, classic);

  }
  
  virtual void TearDown() {
    words_as_paths.clear();
    this->positions_to_read_again.clear();
    this->pvaluesFw.clear();
    this->pvaluesRev.clear();
    this->scoresFw.clear();
    this->scoresRev.clear();
    delete sequences;
    delete matrix;
  }
  
  void pvalues(optimization_type type) {
      std::thread fw_thread = matrix->getPValues(scoresFw, pvaluesFw);
      fw_thread.join();
      std::thread rev_thread = matrix->getPValues(scoresRev, pvaluesRev);
      rev_thread.join();
  }
  
  
  std::vector<char> temp;
  std::vector<double> pvaluesFw, pvaluesRev, scoresFw, scoresRev;
  std::stringstream stream;
  std::set<size_t> positions_to_read_again;
  Pwm* matrix;
  ChromoVector* sequences;
  std::vector<std::vector<char> > words_as_paths;
};





/*
BENCHMARK(InitModel, Baseline, 2, 1) {
  Pwm matrix(global_params.matrix_filename, global_params.p_value);
}


BENCHMARK_F(GetWordsTest, Baseline, 2, 1) {
  this->getWords(classic);
}

BENCHMARK_F(GetWordsTest, TwoLetterAlphabetOptimization, 2, 1) {
  this->getWords(two_letter);
}

BENCHMARK_F(GetWordsTest, FourLetterAlphabetOptimization, 2, 1) {
  this->getWords(four_letter);
}

BENCHMARK_F(InitAhoCAutomata, Baseline, 2, 1) {
  this->initAutomata();
}

BENCHMARK_F(SearchAhoC, SinglePortionBaseline, 2, 1) {
  for (size_t s_i = 0; s_i < this->sequences->size(); s_i++) {
    for (size_t p_i = 0; p_i < this->sequences->getNumberOfParts(s_i); p_i++) {
      ahoc_search(s_i, 
                  p_i,
                  *(this->atm),
                  *(this->matrix),
                  *(this->sequences), 
                  std::string("/tmp/"),
                  std::string("bsf-bench"),
                  this->total_words,
                  *(this->files_to_merge));
    }
  }
}

BENCHMARK_F(RecalcScores, Baseline, 2, 1) {
  this->recalc_score(classic);
}

BENCHMARK_F(RecalcScores, TwoLetterAlphabetOptimization, 2, 1) {
  this->recalc_score(two_letter);
}

BENCHMARK_F(RecalcScores, FourLetterAlphabetOptimization, 2, 1) {
  this->recalc_score(four_letter);
}

*/

/*
BENCHMARK_F(ReadWords, Baseline, 2, 1) {
  this->read_old();
}
*/
/*
BENCHMARK_F(ReadWords, VMM, 1, 1) {
  this->read_new();
}
*/

BENCHMARK_F(CalcPvalues, Baseline, 1, 1) {
  this->pvalues(classic);
}

BENCHMARK_F(CalcPvalues, Advanced, 1, 1) {
  this->pvalues(distance);
}

int main() {
  hayai::ConsoleOutputter consoleOutputter;

  hayai::Benchmarker::AddOutputter(consoleOutputter);
  
  
  hayai::Benchmarker::RunAllTests();
  
  
  
//  return predict(mem_allowed, matrix_filename, p_value, result_folder, result_filename, status_folder, status_filename, filenames, chromonames, starts, ends);
  return 0;
}
