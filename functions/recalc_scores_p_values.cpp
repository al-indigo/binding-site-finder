#include <set>
#include <iostream>
#include <istream>
#include <cstdio>
#include "../structs/chromovector.h"
#include "../structs/pwm.h"

void recalc_scores_p_values(std::istream& fin, 
                            std::vector<size_t>& positions_to_read_again, 
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
    positions_to_read_again.push_back(buf);
  }
  if (positions_to_read_again.empty()) {
    return;
  }
  std::vector<std::vector<char> > words_as_paths(positions_to_read_again.size(), std::vector<char>(matrix.getLength()));

  size_t counter = 0;
  for (std::vector<size_t>::iterator i = positions_to_read_again.begin(); i != positions_to_read_again.end(); ++i, counter++) {
    std::vector<char> temp(matrix.getLength());
    sequences.getWordAsPathTest(sequence_id, *i, matrix.getLength(), temp);
    words_as_paths[counter] = temp;
  }

  matrix.getScores(words_as_paths, scoresFw, scoresRev);

  pvaluesFw = matrix.getPValues(scoresFw);
  pvaluesRev = matrix.getPValues(scoresRev);

}

void recalc_scores_p_values(FILE * fin, 
                            std::vector<size_t>& positions_to_read_again, 
                            size_t mem_allowed, 
                            Pwm& matrix, 
                            ChromoVector& sequences, 
                            size_t sequence_id, 
                            std::vector<double>& pvaluesFw, 
                            std::vector<double>& pvaluesRev, 
                            std::vector<double>& scoresFw, 
                            std::vector<double>& scoresRev) {  
  size_t buf = 0;
  for (int k = 0; k < mem_allowed / 8 && !feof(fin) && !ferror(fin) && fscanf(fin,"%lx\n", &buf); k++) {
    positions_to_read_again.push_back(buf);
  }
  if (positions_to_read_again.empty()) {
    return;
  }
/*  std::vector<std::vector<char> > words_as_paths(positions_to_read_again.size(), std::vector<char>(matrix.getLength()));

  size_t counter = 0;
  for (std::vector<size_t>::iterator i = positions_to_read_again.begin(); i != positions_to_read_again.end(); ++i, counter++) {
    std::vector<char> temp(matrix.getLength());
    sequences.getWordAsPathTest(sequence_id, *i, matrix.getLength(), temp);
    words_as_paths[counter] = temp;
  } */

//  matrix.getScores(words_as_paths, scoresFw, scoresRev);
  std::pair<double, double> scores;
  for (std::vector<size_t>::iterator i = positions_to_read_again.begin(); i != positions_to_read_again.end(); ++i) {
    if (sequences.getWordScores(sequence_id, *i, matrix, scores)) {
      scoresFw.push_back(scores.first);
      scoresRev.push_back(scores.second);
    }
  }

  pvaluesFw = matrix.getPValues(scoresFw);
  pvaluesRev = matrix.getPValues(scoresRev);

}
