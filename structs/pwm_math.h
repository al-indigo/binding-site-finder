#ifndef PWM_MATH_H
#define PWM_MATH_H

#include "pwm_matrix.h"

#include <map>
#include <vector>

#include <limits>
#include <cmath>


#define DISCRETIZATION_VALUE 10000

//NOTE: This is just ported version from perfectosape; no optimizations

//TODO: think about algorithm to make it bulk operational
//TODO: NOTE: no support for backgrounds for now (A=25%, C=25%, G=25%, T=25%)

/*
class own_double_less : public std::binary_function<double,double,bool>
{
public:
  bool operator()( const double &left, const double &right  ) const
  {
    return (fabs(left - right) > std::numeric_limits<double>::epsilon()*DISCRETIZATION_VALUE && left < right);
  }
};
*/
#define own_double_less std::less<double>

double threshold_by_pvalue (double p_value, double * scores_optimistic, pwmMatrix& matrix, std::map<double, double, own_double_less>& distribution);

void count_distribution_after_threshold (pwmMatrix& matrix, double * scores_optimistic, double threshold, std::map<double, double, own_double_less>& distribution);

void pvalues_by_thresholds(std::vector<double>& thresholds, double cutoff, double * scores_optimistic, pwmMatrix& matrix, std::vector<double>& p_values, std::map<double, double, own_double_less>& distribution);


#endif