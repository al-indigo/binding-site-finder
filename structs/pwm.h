#ifndef PWM_H
#define PWM_H

#include "pwm_matrix.h"
#include "pwm_path.h"
#include "pwm_math.h"

#include <cstring>
#include <string>


class Pwm {
  std::string pwmName;
  pwmMatrix pwmFw;
  pwmMatrix pwmRev;
  pwmMatrix pwmCeiled;
//  pwmMatrix pwmFloored;
  
  double threshold;
  
  // We need this scores to understand if current path could contain interesting values
  double * optimisticScoresFw;
  double * optimisticScoresRev;
  double * optimisticScoresCeiledFw;
//  double * optimisticScoresFlooredFw;
  
  // This scores are needed to understand if we have too much variants. For the future
//  double * pessimisticScoresFw;
//  double * pessimisticScoresRev;
  
  pwmPath lastPath;
  
public:
    Pwm (std::string pwmFilename, double p_value);// double p-value
    
    ~Pwm ();
    
    std::vector <double>                threshold_by_pvalues (std::vector<double> pvalues);
    size_t                              getLength() { return lastPath.getLength(); }
    double *                            InitScoresAheadOptimistic(pwmMatrix & pwm);
    std::vector<std::vector<char> >     getWords(unsigned int count);
    bool                                hasMoreWords() { return lastPath.final; };
    std::vector<double>                 getPValues(std::vector<double>& thresholds);
    /* Maybe we will use it to predict the worst case; not now, it's overkill for now */
    //void InitScoresAheadPessimistic(double * scoreVector, pwmMatrix & pwm);
};

#endif
