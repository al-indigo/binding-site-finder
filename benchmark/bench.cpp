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
#include "../utils/jsonxx.h"

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

void parse_options(const char * task_filename, 
                   size_t& mem_allowed, 
                   size_t& genome_len, 
                   std::string& matrix_filename, 
                   double& p_value, 
                   std::set<std::string>& what_to_test
                   ) {
  std::ifstream options(task_filename);
  if (!options) {
    std::cerr << "Task doesn't exist: " << task_filename << std::endl;
    exit(1);
  }
  
  jsonxx::Object task;
  assert(task.parse(options));
  
  mem_allowed = task.get<jsonxx::Number>("memory");
  genome_len = task.get<jsonxx::Number>("genome-len");
  matrix_filename = task.get<jsonxx::String>("matrix");
  p_value = task.get<jsonxx::Number>("p-value");
  
  for (size_t i = 0; i < task.get<jsonxx::Array>("tests").size(); i++) {
    what_to_test.insert(task.get<jsonxx::Array>("tests").get<jsonxx::String>(i));
  }
}

std::string get_working_path()
{
   char path[1024];   
#ifdef __APPLE__
   uint32_t size = sizeof(path);
   if (_NSGetExecutablePath(path, &size) == 0) {
     return std::string(path);
   } else {
     exit(-1);
   }
#endif
   
#ifdef __linux
   return std::string("linux");
#endif
   
   return std::string("nothing");
   
}


int main() {
  size_t mem_allowed;
  size_t genome_len;
  std::string matrix_filename;
  double p_value;
  std::string result_folder;
  std::string result_filename;
  std::string status_folder;
  std::string status_filename;
  std::set<std::string> what_to_test;
  
  std::vector<std::string> filenames;
  std::vector<std::string> chromonames;
  std::vector<size_t> starts;
  std::vector<size_t> ends;
  

  parse_options("benchmark.json", mem_allowed, genome_len, matrix_filename, p_value, what_to_test);
  
  std::cout << "Testing the following: " << std::endl;
  std::cout << "\tMemory allowed to use: " << mem_allowed << std::endl;
  std::cout << "\tModel is: " << matrix_filename << std::endl;
  std::cout << "\tSelected p_value: " << p_value << std::endl;
  std::cout << "\tGenome lenght to generate: " << genome_len << std::endl;
  std::cout << "\tTasks to test: " << std::endl;
  
  for (std::set<std::string>::iterator i = what_to_test.begin(); i != what_to_test.end(); ++i) {
    std::cout << "\t\t" << *i << std::endl;
  }
  
  if (what_to_test.count(std::string("compare"))) {
    std::cout << "\tCompare task is present so making each test for each implementation" << std::endl;
  }
  
  std::string chromofile = get_working_path() + std::string("-chr_random.plain");
  std::ofstream chromo(chromofile);
  static const char alphabet[] = "ACGT";
  for (size_t i=0; i < genome_len; i++) {
    chromo << alphabet[rand() % (sizeof(alphabet)-1)];
  }
  
  chromo.close();
  

  
  
  

  
//  return predict(mem_allowed, matrix_filename, p_value, result_folder, result_filename, status_folder, status_filename, filenames, chromonames, starts, ends);
  return 0;
}