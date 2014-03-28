#ifndef PWM_H
#define PWM_H

#include "pwm_matrix.h"
#include "pwm_path.h"
#include "pwm_math.h"

#include <cstring>
#include <string>
#include <utility>
#include <stdint.h>

enum optimization_type {classic, two_letter, four_letter, distance};

class Pwm {
  std::string pwmName;
  pwmMatrix pwmFw;
  pwmMatrix pwmRev;
  pwmMatrix pwmCeiled;
  
  double threshold;
  
  // We need this scores to understand if current path could contain interesting values
  double * optimisticScoresFw;
  double * optimisticScoresRev;
  double * optimisticScoresCeiledFw;

  std::vector<std::map<uint32_t, double> > precalcmapFw4letter;
  std::vector<std::map<uint32_t, double> > precalcmapRev4letter;
  std::vector<std::map<uint16_t, double> > precalcmapFw2letter;
  std::vector<std::map<uint16_t, double> > precalcmapRev2letter;
  
  pwmPath lastPath;
  
private:
    void                                initPrecalcMap (std::vector<std::map<uint32_t, double> >& precalcmap, pwmMatrix& pwm);
    void                                initPrecalcMap (std::vector<std::map<uint16_t, double> >& precalcmap, pwmMatrix& pwm);
    
    std::vector<std::vector<char> >     getWordsClassic(unsigned int count);
    std::vector<std::vector<char> >     getWords2letter(unsigned int count);
    std::vector<std::vector<char> >     getWords4letter(unsigned int count);
    
    void                                getScoresClassic(std::vector<std::vector<char> >& words, std::vector<double>& scoresFw, std::vector<double>& scoresRev);
    void                                getScores2letter(std::vector<std::vector<char> >& words, std::vector<double>& scoresFw, std::vector<double>& scoresRev);
    void                                getScores4letter(std::vector<std::vector<char> >& words, std::vector<double>& scoresFw, std::vector<double>& scoresRev);
  
public:
    Pwm (std::string pwmFilename, double p_value);            // pwm_type for optimisations via 2-byte alphabet and 4-byte alphabet.
    
    ~Pwm ();
    
    double                              get_threshold_by_pvalue (double pvalue);
    size_t                              getLength() { return lastPath.getLength(); }
    double *                            InitScoresAheadOptimistic(pwmMatrix & pwm);
    std::vector<std::vector<char> >     getWords(unsigned int count, optimization_type _type = classic);

    bool                                hasMoreWords() { return lastPath.final; };
    std::vector<double>                 getPValues(std::vector<double>& thresholds, optimization_type _type = distance);
    
    void                                getScores(std::vector<std::vector<char> >& words, std::vector<double>& scoresFw, std::vector<double>& scoresRev, optimization_type _type = classic);


};

#endif
