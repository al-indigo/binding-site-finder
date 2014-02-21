#ifndef PWM_H
#define PWM_H

#include "pwm_matrix.h"
#include "pwm_path.h"

#include <cstring>
#include <string>


class Pwm {
  std::string pwmName;
  pwmMatrix pwmFw;
  pwmMatrix pwmRev;
  
  double threshold;
  
  // We need this scores to understand if current path could contain interesting values
  volatile double * optimisticScoresFw;
  volatile double * optimisticScoresRev;
  
  // This scores are needed to understand if we have too much variants. For the future
//  double * pessimisticScoresFw;
//  double * pessimisticScoresRev;
  
  pwmPath lastPath;
  
public:
    Pwm (std::string pwmFilename, double threshold);
    ~Pwm ();
    
    size_t                              getLength() { return lastPath.getLength(); }
    double *                            InitScoresAheadOptimistic(pwmMatrix & pwm);
    std::vector<std::vector<char> >     getWords(double threshold, unsigned int count);
    bool                                hasMoreWords() { return lastPath.final; };
    /* Maybe we will use it to predict the worst case; not now, it's overkill for now */
    //void InitScoresAheadPessimistic(double * scoreVector, pwmMatrix & pwm);
};

#endif
