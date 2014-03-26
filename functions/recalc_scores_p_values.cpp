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
  for (int k = 0; k < mem_allowed / 8 && fin >> buf; k++) {
    positions_to_read_again.insert(buf);
  }
  if (positions_to_read_again.empty()) {
    return;
  }
  std::vector<std::vector<char> > words_as_paths(positions_to_read_again.size(), std::vector<char>(matrix.getLength()));

  size_t counter = 0;
  for (std::set<size_t>::iterator i = positions_to_read_again.begin(); i != positions_to_read_again.end(); ++i, counter++) {
    std::vector<char> temp(matrix.getLength());
    sequences.getWordAsPathTest(sequence_id, *i, matrix.getLength(), temp);
    words_as_paths[counter] = temp;
  }

  matrix.getScores(words_as_paths, scoresFw, scoresRev);

  pvaluesFw = matrix.getPValues(scoresFw);
  pvaluesRev = matrix.getPValues(scoresRev);

}