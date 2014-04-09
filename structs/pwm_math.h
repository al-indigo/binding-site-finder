#ifndef PWM_MATH_H
#define PWM_MATH_H

#include "pwm_matrix.h"

#include <map>
#include <vector>

#define DISCRETIZATION_VALUE 10000

//NOTE: This is just ported version from perfectosape; no optimizations

//TODO: think about algorithm to make it bulk operational
//TODO: NOTE: no support for backgrounds for now (A=25%, C=25%, G=25%, T=25%)

double threshold_by_pvalue (double p_value, double * scores_optimistic, pwmMatrix& matrix);

void pvalues_by_thresholds(std::vector<double>& thresholds, double cutoff, double * scores_optimistic, pwmMatrix& matrix, std::vector<double>& p_values);
void pvalues_by_thresholds2(std::vector<double>& thresholds, double cutoff, double * scores_optimistic, pwmMatrix& matrix, std::vector<double>& p_values);


#endif