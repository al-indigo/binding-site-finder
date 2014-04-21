#ifndef PWM_H
#define PWM_H

#include "pwm_matrix.h"
#include "pwm_path.h"
#include "pwm_math.h"

#include <cstring>
#include <string>
#include <utility>
#include <stdint.h>
#include <thread>

enum optimization_type {classic, two_letter, four_letter, distance};
//TODO: make matrix type clearer.

class Pwm {
  std::string pwmName;
  pwmMatrix pwmFw;
  pwmMatrix pwmRev;
  pwmMatrix pwmCeiled;
  
  double threshold;
  size_t words_to_find;
  size_t words_found;
  
  // We need this scores to understand if current path could contain interesting values
  double * optimisticScoresFw;
  double * optimisticScoresRev;
  double * optimisticScoresCeiledFw;

  pwmPath lastPath;
  
  matrix_optimization_type type;
  
  std::map<double, double> distribution_after_threshold;
  
private:
    void                                initPrecalcMap (std::vector<std::map<uint32_t, double> >& precalcmap, pwmMatrix& pwm);
    void                                initPrecalcMap (std::vector<std::map<uint16_t, double> >& precalcmap, pwmMatrix& pwm);
    
    std::vector<std::vector<char> >     getWordsClassic(unsigned int count);

    void                                fillMatrixPwmClassic(std::vector<double>& scorebuffer);
    void                                fillMatrixPwmBit(std::vector<double>& scorebuffer);
//    void                                fillMatrixMoods();
  
public:
    Pwm (std::string pwmFilename, double p_value, matrix_optimization_type _type=pat);
   ~Pwm ();
    
    double                              get_threshold_by_pvalue (double pvalue);
    size_t                              getLength() { return lastPath.getLength(); }
    double *                            InitScoresAheadOptimistic(pwmMatrix & pwm);
    std::vector<std::vector<char> >     getWords(unsigned int count, optimization_type _type = classic);

    bool                                hasMoreWords() { return !lastPath.final; };
    
    std::thread                         getPValuesThreaded(std::vector<double>& thresholds, std::vector<double>& p_values, optimization_type _type = distance);
    void                                getPValuesPlain(std::vector<double>& thresholds, std::vector<double>& p_values, optimization_type _type = distance);
    
    void                                getScores(std::vector<std::vector<char> >& words, std::vector<double>& scoresFw, std::vector<double>& scoresRev, optimization_type _type = classic);
    
    double                        getThreshold() {return this->threshold;}
    size_t                        getNumberOfWords() {return this->words_to_find; };
    size_t                        getNumberOfWordsFound() {return this->words_found; };
    unsigned int                  getNumberOfWordLeft() {return this->words_to_find - this->words_found; };
    matrix_optimization_type      getMatrixType() {return this->type; }
    
    bool                          getScorePair (char * word, std::pair<double, double>& scores);
    void                          getScoresVector (char* buffer, std::vector< bool >& need_to_check, std::vector< double >& scores, std::vector< bool >& strand, std::vector< uint32_t >& positions, size_t offset);


};

#endif
