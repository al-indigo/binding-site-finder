#ifndef PWM_MATRIX_H
#define PWM_MATRIX_H

typedef struct pwm_matrix {
  double * matrix;
  size_t rows;
  size_t cols;
  pwm_matrix() { matrix = NULL; }
 ~pwm_matrix() { delete [] matrix; }
  void init(size_t width, size_t height) { matrix = new double[width*height]; cols = width; rows = height; }
  double& operator()(int row, int col) { return matrix[col + row*cols]; }
} pwmMatrix;

#endif