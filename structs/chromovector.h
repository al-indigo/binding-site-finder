#ifndef CHROMOVECTOR_H
#define CHROMOVECTOR_H

#include "chromo.h"

#include <vector>

class ChromoVector {
  std::vector<Chromo*> chromovector;
  std::vector<std::vector<std::string> > tempresults;
  std::vector<std::string> result_filenames;

public:
  ChromoVector (std::vector<std::string> filenames,
                std::vector<std::string> description,
                std::vector<size_t> start,
                std::vector<size_t> end,
                size_t intersection_size
               );
  size_t size();
  std::string& getFilename(size_t sequence_number);
  std::string& getDescription(size_t sequence_number);
  char * getSeq (size_t sequence_number, size_t part_number);
  size_t getNumberOfParts (size_t sequence_number);
  size_t getPartLength ( size_t sequence_number, size_t part_number );
  size_t getSeqLength (size_t sequence_number);
  size_t getAbsoluteOffset (size_t sequence_number, size_t part_number);
  void releasePart ( size_t sequence_number, size_t part_number);
  
  void getWordsAsPaths (size_t sequence_number, std::vector<size_t>& positions, size_t length, std::vector <std::vector<char> >& result) { chromovector[sequence_number]->getWordsAsPaths(positions, length, result); };
  void getWordAsPathTest (size_t sequence_number, size_t position, size_t length, std::vector<char>& result) { chromovector[sequence_number]->getWordAsPathTest(position, length, result); }
  bool getWordScores (size_t sequence_number, size_t position, Pwm& matrix, std::pair<double, double>& scores) { return chromovector[sequence_number]->getWordScores(position, matrix, scores); }
  void getManyWordScores (size_t sequence_number, std::vector<size_t>& positions_to_read_again, Pwm& matrix, std::vector<double>& scores, std::vector<size_t>& matched, std::vector<bool>& strand) { return chromovector[sequence_number]->getManyWordScores(positions_to_read_again, matrix, scores, matched, strand); }
  void getWordScoresVector (size_t sequence_number, size_t first, size_t last, Pwm& matrix, std::vector<double>& scores, std::vector<bool>& strand, std::vector<size_t>& positions) { return chromovector[sequence_number]->getWordScoresVector(first, last, matrix, scores, strand, positions); }

 ~ChromoVector ( );
};

#endif // CHROMOVECTOR_H
