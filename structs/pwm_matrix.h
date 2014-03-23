#ifndef PWM_MATRIX_H
#define PWM_MATRIX_H

#include <limits>
#include <cstdlib>
#include <string.h>

typedef struct pwm_matrix {
  double * matrix;
  size_t rows;
  size_t cols;
  pwm_matrix() { matrix = NULL; }
 ~pwm_matrix() { delete [] matrix; }
  
  /* NOTE: we are initializing rows to be % 4. Not needed will be 0, doesn't matter.
   * NOTE: DANGER: I know that initing double array with 0 will give us some inaccuracy, but it should be pretty small to notice.
   */
  void init(size_t width, size_t height) { 
//    matrix = new double[width*height];
    matrix = new double[width*height*((4 + height -1)/4)]; 
    cols = width; 
    rows = height; 
    memset(matrix, 0, sizeof(double) * width*height*((4 + height -1)/4)); 
  }  
  
  double& operator()(int row, int col) { return matrix[col + row*cols]; }
} pwmMatrix;

#endif
