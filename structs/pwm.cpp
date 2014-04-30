#include <limits>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <functional>

#include "pwm.h"

void Pwm::fillMatrixPwmClassic(std::vector<double>& scorebuffer) {
  auto count = scorebuffer.size();
  for (auto k=0; k<count; k++) {
    pwmFw.matrix[k] = scorebuffer[k]; //test[k];
    pwmRev.matrix[k] = scorebuffer[count-k-1]; //test[count-k-1];
    pwmCeiled.matrix[k] = ceil(scorebuffer[k] * DISCRETIZATION_VALUE);
  }
}


void Pwm::fillMatrixPwmBit(std::vector<double>& scorebuffer) {
  auto count = scorebuffer.size();

    for (auto i=0; i < count/4; i++) {
      pwmFw.matrix[i*4] = scorebuffer[i*4];
      pwmFw.matrix[i*4+1] = scorebuffer[i*4+1];
      pwmFw.matrix[i*4+2] = scorebuffer[i*4+3];
      pwmFw.matrix[i*4+3] = scorebuffer[i*4+2];
      
      pwmRev.matrix[i*4] = scorebuffer[count - i*4-1];
      pwmRev.matrix[i*4+1] = scorebuffer[count - i*4 - 2];
      pwmRev.matrix[i*4+2] = scorebuffer[count - i*4 - 4];
      pwmRev.matrix[i*4+3] = scorebuffer[count - i*4 - 3];
    } 
    for (auto k=0; k<count; k++) {
      pwmCeiled.matrix[k] = ceil(pwmFw.matrix[k] * DISCRETIZATION_VALUE);
    }
}


Pwm::Pwm(std::string pwmFilename, double p_value, matrix_optimization_type _type): words_found(0), type(_type) {
  std::ifstream ifs(pwmFilename.c_str());
  std::string tempStr;
  size_t count = 0;
  std::vector<double> scorebuffer;
  double dummy;

  if (!(ifs >> std::noskipws >> this->pwmName)) {
    std::cout << "Wrong input: matrix file should have a name at first string" << std::endl;
    exit(-1);
  }
#ifdef DDEBUG_PRINT
  std::cout << "Read matrix name: " << this->pwmName << std::endl;
#endif

  while(ifs >> std::skipws >> dummy) {
    if (!std::isnan(dummy)) {
      scorebuffer.push_back(dummy);
//      test[count] = dummy;
      count++;
    } else {
      std::cout << "Wrong input: matrix can not be parsed as doubles" << std::endl;
    }

  }
#ifdef DDEBUG_PRINT
  std::cout << "Have read " << count << " doubles as matrix content" << std::endl;
#endif

  pwmFw.init(NUM_COLS, count/NUM_COLS);
  pwmRev.init(NUM_COLS, count/NUM_COLS);
  pwmCeiled.init(NUM_COLS, count/NUM_COLS);
//  pwmFloored.init(NUM_COLS, count/NUM_COLS);
  switch (type) {
    case pat: 
              fillMatrixPwmClassic(scorebuffer);
              break;
    case pat_bit_optimization:
              fillMatrixPwmBit(scorebuffer);
              break;
  }
  pwmFw.fillq();
  pwmRev.fillq();
//  std::thread fw(&pwmMatrix::fillq, std::ref(pwmFw));
//  std::thread rev(&pwmMatrix::fillq, std::ref(pwmRev));

#ifdef DDEBUG_PRINT
  std::cout << "Fw Q matrix is " << std::endl;
  pwmFw.printqMatrix();
  std::cout << "Rev Q matrix is " << std::endl;
  pwmRev.printqMatrix();
#endif
  
  lastPath.init(pwmFw.rows);

#ifdef DDEBUG_PRINT
  std::cout.precision(17);
  std::cout << "Forward matrix is:" << std::endl;
  for (int n=0; n < pwmFw.rows; n++) {
    for (int m=0; m < pwmFw.cols; m++) {
      std::cout << pwmFw(n,m) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "Reverse matrix is:" << std::endl; 
  for (int n=0; n < pwmRev.rows; n++) {
    for (int m=0; m < pwmRev.cols; m++) {
      std::cout << pwmRev(n,m) << " ";
    }
    std::cout << std::endl;
  }
#endif

  optimisticScoresFw = InitScoresAheadOptimistic(pwmFw);
  optimisticScoresRev = InitScoresAheadOptimistic(pwmRev);
  optimisticScoresCeiledFw = InitScoresAheadOptimistic(pwmCeiled);
  
  threshold = this->get_threshold_by_pvalue(p_value);
  words_to_find =  2 * (size_t)(floor(p_value * pow(4, getLength())) - 1);
  
//  fw.join();
//  rev.join();

  return;
}


double Pwm::get_threshold_by_pvalue ( double pvalue ) {
  return threshold_by_pvalue(pvalue, optimisticScoresCeiledFw, pwmCeiled, distribution_after_threshold);
}


std::thread Pwm::getPValuesThreaded(std::vector<double>& thresholds, std::vector<double>& p_values, optimization_type _type) {
  return std::thread(pvalues_by_thresholds, std::ref(thresholds), std::ref(threshold), std::ref(optimisticScoresCeiledFw), std::ref(pwmCeiled), std::ref(p_values), std::ref(distribution_after_threshold));
}

void Pwm::getPValuesPlain(std::vector<double>& thresholds, std::vector<double>& p_values, optimization_type _type) {
  pvalues_by_thresholds(thresholds, threshold, optimisticScoresCeiledFw, pwmCeiled, p_values, distribution_after_threshold);
}


double * Pwm::InitScoresAheadOptimistic(pwmMatrix & pwm) {
  double * scoreVector = new double[pwm.rows + 1];
  memset(scoreVector, 0,  sizeof(double) * (pwm.rows + 1));
  for (int i = pwm.rows - 1; i > 0; --i) {
    double max = -1 * std::numeric_limits< double >::max();
    for (unsigned int j = 0; j < pwm.cols; ++j) {
      max = max >= pwm(i,j) ? max : pwm(i,j);
    }
    scoreVector[i-1] = scoreVector[i] + max;
  }

#ifdef DDEBUG_PRINT
  for (int i = 0; i < pwm.rows; i++) {
    std::cout << scoreVector[i] << " ; " << std::endl;
  }
#endif

  return scoreVector;
}


std::vector<std::vector<char> > Pwm::getWords(unsigned int count, optimization_type type) {
  switch (type) {
    case classic:               return getWordsClassic(count);
    case two_letter:            
    case four_letter:           
    default:                    std::cerr << "Smth really wrong: no optimization type at all or optimization is not implemented" << std::endl; exit(-10);
  }
}


std::vector<std::vector<char> > Pwm::getWordsClassic(unsigned int count) {
  
  std::vector<std::vector<char> > words;
  words.resize(std::min(count, count), std::vector<char>(this->getLength()+1));
//  words.resize(std::min(count, this->getNumberOfWordLeft()), std::vector<char>(this->getLength()+1));
  
  std::cout << " words vector size is " << std::min(count, this->getNumberOfWordLeft()) << " , patterns allowed: " << count << std::endl;
  
  /* wordcount is not a regular counter: it follows fit words, not number of iterations */
  
  size_t fw_fits = 0;
  size_t rev_fits = 0;
  
#ifdef DDEBUG_PRINT
  std::map< std::string, double> check; 
#endif
  
  size_t wordcount = 0;
  for (; wordcount < count && !lastPath.final; lastPath.incr()) {
    double fwScore = 0.0;
    double revScore = 0.0;
    bool needBreak = false;
    for (auto depth = 0; depth < lastPath.length; depth++) {
      fwScore += pwmFw(depth, lastPath.path[depth]);
      revScore += pwmRev(depth, lastPath.path[depth]);
      /* here we check if we could ever get more scores on current path */
      if(fwScore + optimisticScoresFw[depth] < threshold &&
         revScore + optimisticScoresRev[depth] < threshold) {
            needBreak = true;
            /* if we can not; we increase path[depth] by one and null all after that */
            lastPath.makeHop(depth);
            break;
      }
    }
    if (needBreak) continue;

    if (fwScore >= threshold || revScore >= threshold) {
      words[wordcount] = lastPath.getWord(type);
      wordcount++;
      words_found++;
#ifdef DDEBUG_PRINT
      check[std::string(&(lastPath.getWord(type)[0]))] = (fwScore > threshold ? fwScore : revScore);
#endif
      fwScore > threshold ? fw_fits++ : rev_fits++;
    }
  }
  words.resize(wordcount);
  

  std::cout << "Internal check: #patterns fit: " << words.size() << " fw: " << fw_fits << " rev: " << rev_fits << std::endl;
#ifdef DDEBUG_PRINT
  for (std::map<std::string, double>::iterator i = check.begin(); i!=check.end(); ++i) {
    std::cout << (*i).first << "\t" << (*i).second << "\n";
  }
#endif

  return words;
}


bool Pwm::getScorePair (char * word, std::pair<double, double>& scores) {
    scores.first = 0.0;
    scores.second = 0.0;
    if (SUPERALPHABET_SIZE > 1) {
      pwmFw.getScoreQ(word, scores.first);
      pwmRev.getScoreQ(word, scores.second);
    } else {
      for (unsigned int k = 0; k < this->getLength(); k++) {
        scores.first += pwmFw(k, word[k]);
        scores.second += pwmRev(k, word[k]);
      }
    }
      
    return (scores.first > threshold || scores.second > threshold)? true: false;
}

void Pwm::getScores (char * word, std::vector<double>& scores, std::vector<bool>& strand, std::vector<size_t>& positions, size_t position) {
    double fwScore = 0.0;
    double revScore = 0.0;
    if (SUPERALPHABET_SIZE > 1) {
      pwmFw.getScoreQ(word, fwScore);
      pwmRev.getScoreQ(word, revScore);
    } else {
      for (auto k = 0; k < this->getLength(); k++) {
        fwScore += pwmFw(k, word[k]);
        revScore += pwmRev(k, word[k]);
      }        
    }
    if (fwScore > threshold) {
      scores.push_back(fwScore);
      positions.push_back(position);
      strand.push_back(true);
    }
    if (revScore > threshold) {
      scores.push_back(revScore);
      positions.push_back(position);
      strand.push_back(false);
    }
}


void Pwm::getScoresVector (char * buffer, std::vector<bool>& need_to_check, std::vector<double>& scores, std::vector<bool>& strand, std::vector<size_t>& positions, size_t offset) {
  for (auto i = 0; i < need_to_check.size(); i++) {
    if (need_to_check[i]) {
      char * word = buffer + i;
      double fwScore = 0.0;
      double revScore = 0.0;
      
      if (SUPERALPHABET_SIZE > 1) {
        pwmFw.getScoreQ(word, fwScore);
        pwmRev.getScoreQ(word, revScore);
        
        #ifdef DDEBUG_PRINT      
        double fwScoreTest = 0.0;
        pwmFw.getScoreSimple(word, fwScoreTest);      
        if (abs(fwScore - fwScoreTest) > 0.01) {
          fwScore = 0.0;
          pwmFw.getScoreQ(word, fwScore);
        }
        #endif      
      } else {
        for (auto k = 0; k < this->getLength(); k++) {
          fwScore += pwmFw(k, word[k]);
          revScore += pwmRev(k, word[k]);
        }        
      }
      
      if (fwScore > threshold) {
        scores.push_back(fwScore);
        positions.push_back(i + offset);
        strand.push_back(true);
      }
      if (revScore > threshold) {
        scores.push_back(revScore);
        positions.push_back(i + offset);
        strand.push_back(false);
      }
    }
  }
}


Pwm::~Pwm() {
  delete [] optimisticScoresFw;
  delete [] optimisticScoresRev;
  delete [] optimisticScoresCeiledFw;
}

