#include "functions.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

#include "../utils/prepare_filename.h"
#include "../utils/mergesort_infile.h"

void ahoc_search(size_t sequence_id, 
                 size_t part_id,
                 AhoCorasickPlus& atm,
                 Pwm& matrix,
                 ChromoVector& sequences, 
                 const std::string& result_folder,
                 const std::string& result_filename,
                 size_t& total_words,
                 std::vector<std::vector<std::string> >& files_to_merge) {
  
  std::vector<size_t> result_vector;
  result_vector.reserve(READ_BLOCK_SIZE);

  char * seq = sequences.getSeq(sequence_id, part_id);
  size_t length = sequences.getPartLength(sequence_id, part_id);
  size_t offset = sequences.getAbsoluteOffset(sequence_id, part_id);

  AhoCorasickPlus::Match aMatch;
  size_t occurances = 0;
  std::cout << "Results for " << sequences.getFilename(sequence_id) << " as " << sequences.getDescription(sequence_id) << " part #" << part_id << "\n";

  atm.search(seq, length, false);
  while (atm.findNext(aMatch))
  {
    result_vector.push_back(aMatch.position + offset - matrix.getLength());
    occurances++;
  }

  std::cout << "Occurances found: " << occurances << std::endl;
  if (occurances > 0) {
    std::string tempfilename = prepare_filename(result_folder + result_filename, "bsf-part-", (sequence_id+1)*10000000 + part_id);
    std::ofstream fout(tempfilename.c_str());
    char * fbuf = new char[FSTREAM_BUF_SIZE];
    fout.rdbuf()->pubsetbuf(fbuf, FSTREAM_BUF_SIZE);
    fout << std::setbase(16);

    files_to_merge[sequence_id].push_back(tempfilename);
    
    for (std::vector<size_t>::iterator z = result_vector.begin(); z != result_vector.end(); ++z) {
      fout << *z << "\n";
    }
    
    fout.close();
    delete [] fbuf;
  }

  sequences.releasePart(sequence_id, part_id);
}