#include "functions.h"

#include <fstream>

#include "../utils/mergesort_infile.h"
#include "../utils/prepare_filename.h"

void merge_files(int id,
                 const std::string& result_folder,
                 const std::string& result_filename,
                 std::vector<std::vector<std::string> >& files_to_merge,
                 std::vector<std::string>& merged_files) {
    
  if (files_to_merge[id].size() == 0) return;
    
  time_t t = time(NULL);
  std::string tempfilename = prepare_filename(result_folder + result_filename, std::string("-temp-merge-stl-"), 10000000*(id+1) + 0, &t);
  std::ofstream fout(tempfilename.c_str());
  fout.close();
  
  for (unsigned int j = 0; j < files_to_merge[id].size(); j++) {
    std::string read_file = prepare_filename(result_folder + result_filename, std::string("-temp-merge-stl-"), 10000000*(id+1) + j, &t);
    std::string write_file = prepare_filename(result_folder + result_filename, std::string("-temp-merge-stl-"), 10000000*(id+1) + j+1, &t);
    
    merge_sort(files_to_merge[id][j], read_file, write_file);
  }
  merged_files.push_back(prepare_filename(result_folder + result_filename, std::string("-temp-merge-stl-"), 10000000*(id+1) + files_to_merge[id].size(), &t));
}