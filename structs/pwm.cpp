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

//#define DDEBUG_PRINT

Pwm::Pwm(std::string pwmFilename, double p_value): words_found(0) {
  std::ifstream ifs(pwmFilename.c_str());
  std::string tempStr;
  size_t count = 0;
  std::vector<double> scorebuffer;
//  double test[36];
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
  for (int k=0; k<count; k++) {
    pwmFw.matrix[k] = scorebuffer[k]; //test[k];
    pwmRev.matrix[k] = scorebuffer[count-k-1]; //test[count-k-1];
    pwmCeiled.matrix[k] = ceil(scorebuffer[k] * DISCRETIZATION_VALUE);
//    pwmFloored.matrix[k] = floor(scorebuffer[k] * DISCRETIZATION_VALUE);
  }

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
  
  threshold = threshold_by_pvalue(p_value, optimisticScoresCeiledFw, pwmCeiled);
  words_to_find = (size_t) (floor(p_value * pow(4, getLength())) + 0.01);
  
  initPrecalcMap(precalcmapFw4letter, pwmFw);
  initPrecalcMap(precalcmapRev4letter, pwmRev);
  initPrecalcMap(precalcmapFw2letter, pwmFw);
  initPrecalcMap(precalcmapRev2letter, pwmRev);
  
  return;
}

double Pwm::get_threshold_by_pvalue ( double pvalue ) {
  return threshold_by_pvalue(pvalue, optimisticScoresCeiledFw, pwmCeiled);
}


std::thread Pwm::getPValues(std::vector<double>& thresholds, std::vector<double>& p_values, optimization_type _type) {
  if (_type == distance) {
    return std::thread(pvalues_by_thresholds2, std::ref(thresholds), std::ref(threshold), std::ref(optimisticScoresCeiledFw), std::ref(pwmCeiled), std::ref(p_values));
  } else {
    return std::thread(pvalues_by_thresholds, std::ref(thresholds), std::ref(threshold), std::ref(optimisticScoresCeiledFw), std::ref(pwmCeiled), std::ref(p_values));
  }
}


double * Pwm::InitScoresAheadOptimistic(pwmMatrix & pwm) {
  double * scoreVector = new double[pwm.rows + 1];
  memset(scoreVector, 0,  sizeof(double) * (pwm.rows + 1));
  for (int i = pwm.rows - 1; i > 0; --i) {
    double max = -1 * std::numeric_limits< double >::max();
    for (int j = 0; j < pwm.cols; ++j) {
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
    case two_letter:            return getWords2letter(count);
    case four_letter:           return getWords4letter(count);
    default:                    std::cerr << "Smth really wrong: no optimization type at all" << std::endl; exit(-10);
  }
}

std::vector<std::vector<char> > Pwm::getWordsClassic(unsigned int count) {
  
  std::vector<std::vector<char> > words;
  words.resize(std::min(count, this->getNumberOfWordLeft()), std::vector<char>(this->getLength()));
  
//  words.reserve(count * sizeof(std::string::value_type) * lastPath.length );
  /* wordcount is not a regular counter: it follows fit words, not number of iterations */
  
  size_t fw_fits = 0;
  size_t rev_fits = 0;
  
  for (int wordcount = 0; wordcount < count && !lastPath.final; lastPath.incr()) {
    double fwScore = 0.0;
    double revScore = 0.0;
    bool needBreak = false;
    for (int depth = 0; depth < lastPath.length; depth++) {
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
//      words.push_back(lastPath.getWord());
      words[wordcount] = lastPath.getWord();
      wordcount++;
      words_found++;
      if (words_found >= words_to_find) lastPath.setFinal();
#ifdef DDEBUG_PRINT
      std::cout << "Word fits: " << lastPath.getWord() << " with score: " << (fwScore > threshold ? fwScore : revScore) << std::endl;
#endif
      fwScore > threshold ? fw_fits++ : rev_fits++;
    }
  }
//#ifdef DDEBUG_PRINT
//  std::cout << "Internal check: #patterns fit: " << words.size() << " fw: " << fw_fits << " rev: " << rev_fits << std::endl;
//#endif
  return words;
}


//NOTE: This implementation fails without failing branches prediction; need to think more.
std::vector<std::vector<char> > Pwm::getWords4letter(unsigned int count) {
  int blocksNumber = (sizeof(uint32_t) + pwmFw.rows - 1)/sizeof(uint32_t);
  std::vector<std::vector<char> > words;
  size_t fw_fits = 0;
  size_t rev_fits = 0;
  for (int wordcount = 0; wordcount < count && !lastPath.final; lastPath.incr()) {
    double fwScore = 0.0;
    double revScore = 0.0;
    bool needBreak = false;    
    for (int i = 0; i < blocksNumber; i++) {
      fwScore += precalcmapFw4letter[i][*((uint32_t * )(lastPath.path + i*sizeof(uint32_t)))];
      revScore += precalcmapRev4letter[i][*((uint32_t * )(lastPath.path + i*sizeof(uint32_t)))];
    }
    if (fwScore >= threshold || revScore >= threshold) {
      words.push_back(lastPath.getWord());
      wordcount++;
#ifdef DDEBUG_PRINT
      std::cout << "Word fits: " << lastPath.getWord() << " with score: " << (fwScore > threshold ? fwScore : revScore) << std::endl;
#endif
      fwScore > threshold ? fw_fits++ : rev_fits++;
    }
  }
//#ifdef DDEBUG_PRINT
  std::cout << "Internal check: #patterns fit: " << words.size() << " fw: " << fw_fits << " rev: " << rev_fits << std::endl;
//#endif
  return words;  
}

std::vector<std::vector<char> > Pwm::getWords2letter(unsigned int count) {
  int blocksNumber = (sizeof(uint16_t) + pwmFw.rows - 1)/sizeof(uint16_t);
  std::vector<std::vector<char> > words;
  size_t fw_fits = 0;
  size_t rev_fits = 0;
  for (int wordcount = 0; wordcount < count && !lastPath.final; lastPath.incr()) {
    double fwScore = 0.0;
    double revScore = 0.0;
    bool needBreak = false;    
    for (int i = 0; i < blocksNumber; i++) {
      fwScore += precalcmapFw2letter[i][*((uint16_t * )(lastPath.path + i*sizeof(uint16_t)))];
      revScore += precalcmapRev2letter[i][*((uint16_t * )(lastPath.path + i*sizeof(uint16_t)))];
    }
    if (fwScore >= threshold || revScore >= threshold) {
      words.push_back(lastPath.getWord());
      wordcount++;
#ifdef DDEBUG_PRINT
      std::cout << "Word fits: " << lastPath.getWord() << " with score: " << (fwScore > threshold ? fwScore : revScore) << std::endl;
#endif
      fwScore > threshold ? fw_fits++ : rev_fits++;
    }
  }
//#ifdef DDEBUG_PRINT
  std::cout << "Internal check: #patterns fit: " << words.size() << " fw: " << fw_fits << " rev: " << rev_fits << std::endl;
//#endif
  return words;  
}

// We use this function to minimize operations for calculation of scores
void Pwm::initPrecalcMap (std::vector<std::map<uint32_t, double> >& precalcmap, pwmMatrix& pwm) {
  int blocksNumber = (sizeof(uint32_t) + pwm.rows - 1)/sizeof(uint32_t); //ceiling for int division (looks creepy)
  for (int i = 0; i < blocksNumber; i++) {
    pwmPath tempPath;
    tempPath.init(sizeof(uint32_t));
    std::map<uint32_t, double> curblockmap;
    for ( ; !tempPath.final; tempPath.incr()) {
      // i'm unrolling the loop by hand here to be sure what's going on.
      double score = pwm(i*sizeof(uint32_t), tempPath.path[0]) + pwm(i*sizeof(uint32_t) + 1, tempPath.path[1]) + pwm(i*sizeof(uint32_t) + 2, tempPath.path[2]) + pwm(i*sizeof(uint32_t) + 3, tempPath.path[3]);
      uint32_t key = *((uint32_t * )tempPath.path);
      curblockmap[key] = score;
    }
    precalcmap.push_back(curblockmap);
  }
}

void Pwm::initPrecalcMap (std::vector<std::map<uint16_t, double> >& precalcmap, pwmMatrix& pwm) {
  int blocksNumber = (sizeof(uint16_t) + pwm.rows - 1)/sizeof(uint16_t); //ceiling for int division (looks creepy)
  for (int i = 0; i < blocksNumber; i++) {
    pwmPath tempPath;
    tempPath.init(sizeof(uint16_t));
    std::map<uint16_t, double> curblockmap;
    for ( ; !tempPath.final; tempPath.incr()) {
      // i'm unrolling the loop by hand here to be sure what's going on.
      double score = pwm(i*sizeof(uint16_t), tempPath.path[0]) + pwm(i*sizeof(uint16_t) + 1, tempPath.path[1]);
      uint16_t key = *((uint16_t * )tempPath.path);
      curblockmap[key] = score;
    }
    precalcmap.push_back(curblockmap);
  }
}


Pwm::~Pwm() {
  delete [] optimisticScoresFw;
  delete [] optimisticScoresRev;
  delete [] optimisticScoresCeiledFw;
}


void Pwm::getScores (std::vector< std::vector< char > >& words, std::vector<double>& scoresFw, std::vector<double>& scoresRev, optimization_type type) {
  switch (type) {
    case classic:               return getScoresClassic(words, scoresFw, scoresRev);
    case two_letter:            return getScores2letter(words, scoresFw, scoresRev);
    case four_letter:           return getScores4letter(words, scoresFw, scoresRev);
    default:                    std::cerr << "Smth really wrong: no optimization type at all" << std::endl; exit(-11);
  }
}

bool Pwm::getScorePair (char * word, std::pair<double, double>& scores) {
    double fwScore = 0.0;
    double revScore = 0.0;
    
    for (int k = 0; k < this->getLength(); k++) {
      fwScore += pwmFw(k, word[k]);
      revScore += pwmRev(k, word[k]);
    }
    scores.first = fwScore;
    scores.second = revScore;
    
    return (fwScore > threshold || revScore > threshold)? true: false;
}


void Pwm::getScoresClassic (std::vector< std::vector< char > >& words, std::vector<double>& scoresFw, std::vector<double>& scoresRev) {
  scoresFw.reserve(words.size());
  scoresRev.reserve(words.size());
  std::cout << "Calculating scores again" << std::endl;
  for (std::vector<std::vector<char> >::iterator i = words.begin(); i != words.end(); ++i) {
    std::vector<char> word = *i;
    double fwScore = 0.0;
    double revScore = 0.0;
    
    for (int k = 0; k < word.size(); k++) {
      fwScore += pwmFw(k, word[k]);
      revScore += pwmRev(k, word[k]);
    }
    scoresFw.push_back(fwScore);
    scoresRev.push_back(revScore);
  }
  std::cout << "Scores ready" << std::endl;
  return;

}

void Pwm::getScores4letter (std::vector< std::vector< char > >& words, std::vector<double>& scoresFw, std::vector<double>& scoresRev) {
  scoresFw.reserve(words.size());
  scoresRev.reserve(words.size());
  std::cout << "Calculating scores again (version 2)" << std::endl;
  
  int blocksNumber = (sizeof(uint32_t) + pwmFw.rows - 1)/sizeof(uint32_t);
  char * buffer = new char[pwmFw.rows * (sizeof(uint32_t) + pwmFw.rows - 1)/sizeof(uint32_t) ]; 
  memset(buffer, 0, sizeof(char) * pwmFw.rows * (sizeof(uint32_t) + pwmFw.rows - 1)/sizeof(uint32_t) );

  for (std::vector<std::vector<char> >::iterator i = words.begin(); i != words.end(); ++i) {
    std::vector<char> word = *i;
    memcpy(buffer, &word[0], word.size() * sizeof(char));
    double fwScore = 0.0;
    double revScore = 0.0;
    for (int i = 0; i < blocksNumber; i++) {
      fwScore += precalcmapFw4letter[i][*((uint32_t * )(buffer + i*sizeof(uint32_t)))];
      revScore += precalcmapRev4letter[i][*((uint32_t * )(buffer + i*sizeof(uint32_t)))];
    }
    scoresFw.push_back(fwScore);
    scoresRev.push_back(revScore);
  }
  delete[] buffer;

  std::cout << "Scores ready" << std::endl;
  return;
}

void Pwm::getScores2letter (std::vector< std::vector< char > >& words, std::vector<double>& scoresFw, std::vector<double>& scoresRev) {
  scoresFw.reserve(words.size());
  scoresRev.reserve(words.size());
  std::cout << "Calculating scores again (version 3)" << std::endl;
  
  int blocksNumber = (sizeof(uint16_t) + pwmFw.rows - 1)/sizeof(uint16_t);
  char * buffer = new char[pwmFw.rows * (sizeof(uint16_t) + pwmFw.rows - 1)/sizeof(uint16_t) ]; 
  memset(buffer, 0, sizeof(char) * pwmFw.rows * (sizeof(uint16_t) + pwmFw.rows - 1)/sizeof(uint16_t) );

  for (std::vector<std::vector<char> >::iterator i = words.begin(); i != words.end(); ++i) {
    std::vector<char> word = *i;
    memcpy(buffer, &word[0], word.size() * sizeof(char));
    double fwScore = 0.0;
    double revScore = 0.0;
    for (int i = 0; i < blocksNumber; i++) {
      fwScore += precalcmapFw2letter[i][*((uint16_t * )(buffer + i*sizeof(uint16_t)))];
      revScore += precalcmapRev2letter[i][*((uint16_t * )(buffer + i*sizeof(uint16_t)))];
    }
    scoresFw.push_back(fwScore);
    scoresRev.push_back(revScore);
  }
  delete[] buffer;

  std::cout << "Scores ready" << std::endl;
  return;
}