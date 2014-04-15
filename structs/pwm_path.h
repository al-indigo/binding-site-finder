#ifndef PWM_PATH_H
#define PWM_PATH_H

#include <vector>
#include <cstring>
#include <iostream>
#include <stdint.h>

#define NUM_COLS 4

//enum matrix_type {pcm, moods}; 
enum matrix_optimization_type {pat, pat_bit_optimization}; //bit optimization type is connected with word-to-path procedure in chromo.cpp. There is a bit-optimized version that need a special form of matrix.

/* Inlined for better performance */

typedef struct pwm_path {
  char * path;
  size_t length;
  bool final;
  std::vector<char> word;
  pwm_path() { path = NULL; length = 0; final = false; }
 ~pwm_path() { delete[] path; }
  void init (size_t _length){ 
    path = new char[_length * (sizeof(uint32_t) + _length - 1)/sizeof(uint32_t) ];
    length = _length; 
    memset(path, 0, sizeof(char) * _length * (sizeof(uint32_t) + _length - 1)/sizeof(uint32_t) );
    word.resize(length+1);
  }
  
  /* NOTE: I'm using here vector instead of char[] to avoid leaks in future: it's simplier to maintain */   
  std::vector<char> getWord(matrix_optimization_type type) {
    word[length] = '\0';
    switch (type) {
      case pat: {
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
        break;
      }
      case pat_bit_optimization: {
        for(int i = 0; i < length; i++) {
          switch (path[i]) {
            case 0: 
              word[i] = 'A';
              break;
            case 1:
              word[i] = 'C';
              break;
            case 2:
              word[i] = 'T';
              break;
            case 3:
              word[i] = 'G';
              break;
            default:
              std::cerr << "Something wrong with the matrix: it seems to have more than 4 columns" << std::endl;
              exit(-1);
          }
        }
        break;
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
          setFinal();
          return;
        }
        path[i] = 0;
        continue;
      }
    }
  }
  
  void setFinal() {
    memset(path, NUM_COLS - 1, sizeof(char) * length);
    final = true;
  }

} pwmPath;

#endif
