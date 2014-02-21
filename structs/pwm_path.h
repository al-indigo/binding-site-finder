#ifndef PWM_PATH_H
#define PWM_PATH_H

#include <vector>

#define NUM_COLS 4


/* Inlined for better performance */

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

#endif