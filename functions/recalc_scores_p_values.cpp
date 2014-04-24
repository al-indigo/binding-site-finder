#include <set>
#include <iostream>
#include <istream>
#include <cstdio>
#include <thread>
#include "../structs/chromovector.h"
#include "../structs/pwm.h"


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
  size_t last = -1;
  for (auto k = 0; k < mem_allowed / 8 && !feof(fin) && !ferror(fin) && fscanf(fin,"%lx\n", &buf); k++) {
    if (last != buf) {
      positions_to_read_again.push_back(buf);
      last = buf;
    }
  }
  if (positions_to_read_again.empty()) {
    return;
  }

  std::pair<double, double> scores;
  scoresFw.resize(positions_to_read_again.size());
  scoresRev.resize(positions_to_read_again.size());
  for (size_t i = 0; i < positions_to_read_again.size(); i++) {
     sequences.getWordScores(sequence_id, positions_to_read_again[i], matrix, scores); 
     scoresFw[i] = scores.first;
     scoresRev[i] = scores.second;
  }

  std::thread fw_thread = matrix.getPValuesThreaded(scoresFw, pvaluesFw);
  std::thread rev_thread = matrix.getPValuesThreaded(scoresRev, pvaluesRev);
  fw_thread.join();
  rev_thread.join();

}


void recalc_scores_p_values(FILE * fin, 
                            size_t mem_allowed, 
                            Pwm& matrix, 
                            ChromoVector& sequences, 
                            size_t sequence_id, 
                            std::vector<double>& pvalues, 
                            std::vector<double>& scores,
                            std::vector<bool>& strand,
                            std::vector<size_t>& matched
                           ) {
  std::vector<size_t> positions_to_read_again;
  size_t buf = 0;
  size_t last = -1;
  for (auto k = 0; k < mem_allowed / 8 && !feof(fin) && !ferror(fin) && fscanf(fin,"%lx\n", &buf); k++) {
    if (last != buf) {
      positions_to_read_again.push_back(buf);
      last = buf;
    }
  }
  if (positions_to_read_again.empty()) {
    return;
  }
  matched.reserve(positions_to_read_again.size());
  scores.reserve(positions_to_read_again.size());
  strand.reserve(positions_to_read_again.size());
  
  sequences.getManyWordScores(sequence_id, positions_to_read_again, matrix, scores, matched, strand); 

  matrix.getPValuesPlain(scores, pvalues);

}