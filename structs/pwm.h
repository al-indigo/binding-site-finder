#ifndef PWM_H
#define PWM_H

#include <cstring>
#include <string>
#include <stack>
#include <vector>

#define NUM_COLS 4


typedef struct pwm_matrix {
  double * matrix;
  size_t rows;
  size_t cols;
  pwm_matrix() { matrix = NULL; }
 ~pwm_matrix() { delete matrix; }
  void init(size_t width, size_t height) { matrix = new double[width*height]; cols = width; rows = height; }
  double& operator()(int row, int col) { return matrix[col + row*cols]; }
} pwmMatrix;

typedef struct pwm_path {
  char * path;
  size_t length;
  bool final;
  pwm_path() { path = NULL; length = 0; final = false; }
 ~pwm_path() { delete[] path; }
  void init (size_t _length){ path = new char[_length]; length = _length; memset(path, 0, sizeof(char) * length); }
  
/* NOTE: I'm using here vector instead of char[] to avoid leaks in future: it's simplier to maintain */   
  std::vector<char> getWord() {
    std::vector<char> word(length+1);
    word[length] = '\0';
    for(int i = 0; i < length; i++) {
      switch (path[i]) {
        case 0: 
          word[i] = 'A';
          break;
        case 1:
          word[i] = 'C';
          break;
        case 2:
          word[i] = 'G';
          break;
        case 3:
          word[i] = 'T';
          break;
        default:
          std::cerr << "Something wrong with the matrix: it seems to have more than 4 columns" << std::endl;
          exit(-1);
      }
    }
    return word;
  }
  
  size_t getLength() { 
    return length; 
  }
  
  void incr() {
    for (int i = length-1; i >= 0; i--) {
      if (path[i] < NUM_COLS - 1) {
          path[i]++;
          return;
      } else {
          if (i == 0) {
            final = true;
            return;
          }
          path[i] = 0;
          continue;
      }
    }
    return;
  }
  
  void makeHop(size_t hopPoint) {
    if (hopPoint == length - 1) return;

    for (int i = hopPoint + 1; i < length - 1; i++) {
      path[i] = 0;
    }

    path[length - 1] = -1;

    for (int i = hopPoint; i >= 0; i--) {
      if (path[i] < NUM_COLS - 1) {
        path[i]++;
        return;
      } else {
        /* we have finished */
        if (i == 0 && path[i] == NUM_COLS - 1) {
          memset(path, NUM_COLS - 1, sizeof(char) * length);
          final = true;
          return;
        }
        path[i] = 0;
        continue;
      }
    }
  }

} pwmPath;

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
    
    size_t getLength() { return lastPath.getLength(); }
    double * InitScoresAheadOptimistic(pwmMatrix & pwm);
    std::vector<std::vector<char> > getWords(double threshold, unsigned int count);
    bool hasMoreWords() { return lastPath.final; };
    /* Maybe we will use it to predict the worst case; not now, it's overkill for now */
    //void InitScoresAheadPessimistic(double * scoreVector, pwmMatrix & pwm);
};

#endif
