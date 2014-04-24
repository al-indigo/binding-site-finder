#ifndef PWM_MATRIX_H
#define PWM_MATRIX_H

#include <thread>
#include <limits>
#include <cstdlib>
#include <string.h>
#include <stdint.h>

//IMPORTANT: NOTE: Q stands for superAlphabet size and it should be known at compile-time. MUST be >=1. If it's equal to 1 then simple alphabet is used.
#define SUPERALPHABET_SIZE 5
//#define DDEBUG_PRINT

#ifdef DDEBUG_PRINT
#include <iostream>
#endif

typedef struct pwm_matrix {
  double * matrix;
  double * qmatrix;
  size_t rows;
  size_t cols;
  size_t qrows;
  size_t qcols;
  pwm_matrix() { matrix = NULL; qmatrix = NULL; }
 ~pwm_matrix() { delete [] matrix; delete [] qmatrix; }
  
  /* NOTE: we are initializing rows to be % 4. Not needed rows will be 0, doesn't matter.
   * NOTE: DANGER: I know that initing double array with 0 will give us some inaccuracy, but it should be pretty small to notice.
   */
  void init(size_t width, size_t height) { 
//    matrix = new double[width*height];
    matrix = new double[width*SUPERALPHABET_SIZE*((SUPERALPHABET_SIZE + height -1)/SUPERALPHABET_SIZE)]; 
    cols = width; 
    rows = height; 
    memset(matrix, 0, sizeof(double) * width*SUPERALPHABET_SIZE*((SUPERALPHABET_SIZE + height -1)/SUPERALPHABET_SIZE)); 
    if (SUPERALPHABET_SIZE > 1) {
      initq();
    }
  }
  
  void initq() {
    if (SUPERALPHABET_SIZE > 1) {
      qcols = cols;
      for (unsigned int i = 1; i < SUPERALPHABET_SIZE; i++) qcols *= cols;
      qrows = ((SUPERALPHABET_SIZE + rows - 1)/SUPERALPHABET_SIZE);
      qmatrix = new double[qrows * qcols];
    }
  }
  
  void fillq() {
    if (SUPERALPHABET_SIZE > 1) {
      for (auto row = 0; row < qrows; row++) {
        for (auto column = 0; column < qcols; column++) {
          for (auto i = 0; i < SUPERALPHABET_SIZE; i++) {
            qmatrix[column + row*qcols] += matrix[(0x0000000000000003&(column>>2*(SUPERALPHABET_SIZE-i-1))) + (SUPERALPHABET_SIZE*row + i)*cols];
          }
        }
      }
    }
  }

  void printqMatrix() {
    #ifdef DDEBUG_PRINT  
    for (auto row = 0; row < qrows; row++) {
      for (auto column = 0; column < qcols; column++) {
        std::cout << this->operator()(row, column, true) << "  ";
      }
      std::cout << std::endl;
    }
    #endif
  }

  void getScoreSimple (char * word, double& score) {
    for (auto k = 0; k < rows; k++) {
      score += matrix[word[k] + k*cols];
    }
  }
 
  void getScoreQ (char * word, double& score) {
    size_t idx = 0;
    for (auto k = 0; k < qrows; k++) {
      idx = 0;
      for (auto i = 0; i < SUPERALPHABET_SIZE; i++) {
        idx |= word[k*SUPERALPHABET_SIZE+i]<<(2*(SUPERALPHABET_SIZE-1-i));
      }

      score += qmatrix[idx + k*qcols];
    }
  }
  
  double& operator()(int row, int col) { return matrix[col + row*cols]; }
  double& operator()(int row, int col, bool overload) { return qmatrix[col + row*qcols]; }
} pwmMatrix;

#endif
