#include <set>
#include <iostream>
#include <istream>
#include "../structs/chromovector.h"
#include "../structs/pwm.h"

void recalc_scores_p_values(std::istream& fin, 
                            std::set<size_t>& positions_to_read_again, 
                            size_t mem_allowed, 
                            Pwm& matrix, 
                            ChromoVector& sequences, 
                            size_t sequence_id, 
                            std::vector<double>& pvaluesFw, 
                            std::vector<double>& pvaluesRev, 
                            std::vector<double>& scoresFw, 
                            std::vector<double>& scoresRev) {  
  size_t buf = 0;
  for (int k = 0; k < mem_allowed * 1024 * 1024 / 8 && fin >> buf; k++) {
    positions_to_read_again.insert(buf);
  }
  if (positions_to_read_again.empty()) {
    return;
  }

  std::vector<std::vector<char> > words_as_paths(positions_to_read_again.size());

  sequences.getWordsAsPaths(sequence_id, positions_to_read_again, matrix.getLength(), words_as_paths);

  double tstart, tstop, ttime;
  tstart = (double)clock()/CLOCKS_PER_SEC;

  matrix.getScores(words_as_paths, scoresFw, scoresRev);

  tstop = (double)clock()/CLOCKS_PER_SEC;
  std::cout << "Got scores for " << tstop - tstart << " seconds" << std::endl; 

  tstart = (double)clock()/CLOCKS_PER_SEC;
  pvaluesFw = matrix.getPValues(scoresFw);
  pvaluesRev = matrix.getPValues(scoresRev);
  tstop = (double)clock()/CLOCKS_PER_SEC;
  std::cout << "Got p-values for " << tstop - tstart << " seconds" << std::endl; 

}