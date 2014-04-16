#include "pwm_math.h"

#include <iostream>
#include <cmath>
#include <numeric>
#include <utility>
#include <stdlib.h>
#include <algorithm>

/* TODO: NOTE: I know that using doubles in comparisons (or find) is not semantically right, need to fix sometime. 
 *             It's implemented this way because when using ulong, results are not valid for some tests 
 *             due to accuracy in conversions. */

//Approximated inverse Gauss error function
inline double inverf (double x) {
  //NOTE: don't want to use c++11 stdrt just for this.
  int sign = (x > 0) ? 1 : ((x < 0) ? -1 : 0);
  x = (x > 0) ? x : -x;
  double a = 0.147; //from wiki: this is more accurate value for approximation
  double tmp = (2.0 / (M_PI * a) + (log(1.0 - x * x)) / 2.0);
  double part0 = tmp * tmp;
  double part = -2.0 / (M_PI * a) - log(1.0 - x * x) / 2.0 + sqrt(-1.0 / a * log(1.0 - x * x) + part0);

  return sign * sqrt(part);
}

inline static std::map<double, double> recalc_score_hash (std::map<double, double> scores, pwmMatrix& matrix, size_t row, double least_sufficient) {
  std::map<double, double> new_scores;
  
  for (std::map <double, double>::iterator i = scores.begin(); i != scores.end(); ++i) {
    double score = (*i).first;
    double count = (*i).second;
    for (unsigned int k = 0; k < matrix.cols; k++) {
      double new_score = score + matrix(row, k);
      if (new_score >= least_sufficient) {
        std::map<double, double>::iterator tempit = new_scores.find(new_score);
        if (tempit != new_scores.end()) {
          (*tempit).second += count;
        } else {
          new_scores.insert (std::make_pair(new_score, count));
        }
      }
    }
  }
  return new_scores;
}

void count_distribution_after_threshold (pwmMatrix& matrix, double * scores_optimistic, double threshold, std::map<double, double>& distribution) {
  distribution.insert(std::make_pair(0.0, 1.0));
  for (unsigned int row = 0; row < matrix.rows; row++) {
    distribution = recalc_score_hash(distribution, matrix, row, threshold - scores_optimistic[row]);
  }
}

inline static double score_mean (pwmMatrix& matrix) {
  double result = 0.0;
  for (unsigned int row = 0; row < matrix.rows; row++) {
    double row_sum = 0.0;
    for (unsigned int col = 0; col < matrix.cols; col++) {
      row_sum += matrix(row,col);
    }
    result += row_sum/4.0;
  }
  return result;
}

inline static double score_variance (pwmMatrix& matrix) {
  double result = 0.0;
  for (unsigned int row = 0; row < matrix.rows; row++) {
    double mean_value = 0.0;
    double mean_square = 0.0;
    
    for (unsigned int col = 0; col < matrix.cols; col++) {
      mean_value += matrix(row,col);
      mean_square += matrix(row,col) * matrix(row,col);
    }
    mean_value = mean_value/4.0;
    mean_square = mean_square/4.0;
    
    result += mean_square - mean_value*mean_value;
  }
  return result;
}


inline double static threshold_gauss_estimation (double pvalue, pwmMatrix& matrix) {
  double sigma = sqrt(score_variance(matrix));
  double n_ = inverf(1.0 - 2.0 * pvalue)*sqrt(2.0);
  return score_mean(matrix) + n_ * sigma;
}


inline static size_t vocabularyVolume (pwmMatrix& matrix) {
  return pow(matrix.cols, matrix.rows);
}


inline static std::map<double,double>  count_distribution_under_pvalue (double max_pvalue, double * scores_optimistic, pwmMatrix& matrix) {
    std::map<double, double> cnt_distribution;
    double values_sum = 0.0;
    double look_for_count = max_pvalue * (double)vocabularyVolume(matrix);

    while (values_sum < look_for_count) {
      double approximate_threshold;
      approximate_threshold = threshold_gauss_estimation(max_pvalue, matrix);
      count_distribution_after_threshold(matrix, scores_optimistic, approximate_threshold, cnt_distribution);
      max_pvalue *= 2; // if estimation counted too small amount of words - try to lower threshold estimation by doubling pvalue
      
      values_sum = 0.0;
      for (std::map<double, double>::iterator i = cnt_distribution.begin(); i != cnt_distribution.end(); ++i) {
        values_sum += (*i).second;
      }
    }

    return cnt_distribution;
}


double threshold_by_pvalue (double p_value, double * scores_optimistic, pwmMatrix& matrix, std::map<double, double>& distribution) {
  double threshold = 0.0;
  distribution = count_distribution_under_pvalue(p_value, scores_optimistic, matrix);
  std::vector<double> values, keys, scores_partial_sums;

  scores_partial_sums.reserve(distribution.size());
  for (std::map<double,double>::reverse_iterator i = distribution.rbegin(); i != distribution.rend(); ++i) {
    values.push_back((*i).second);
    keys.push_back((*i).first);
  }
  
  std::partial_sum(values.begin(), values.end(), std::back_inserter(scores_partial_sums));
  
  std::vector<double>::iterator lower = std::lower_bound(scores_partial_sums.begin(), scores_partial_sums.end(), floor(p_value * vocabularyVolume(matrix)));
  size_t l_d = (size_t)abs(std::distance(scores_partial_sums.begin(), lower));

  threshold = (keys[l_d])/DISCRETIZATION_VALUE;
  
  return threshold;
}


void pvalues_by_thresholds(std::vector< double >& thresholds, double cutoff, double* scores_optimistic, pwmMatrix& matrix, std::vector<double>& p_values, std::map<double, double>& distribution) {
  double volume = (double)vocabularyVolume(matrix);
  
  //NOTE: for p_values, values and keys full space allocated at creation. We can't do this for counts_partial_sums becaus later use back_inserter
  std::vector<double> values(distribution.size()), 
                      keys(distribution.size()), 
                      counts_partial_sums;
  p_values.resize(thresholds.size());
  values.resize(distribution.size());
  keys.resize(distribution.size());
                      
  counts_partial_sums.reserve(distribution.size());
  
  size_t aux_counter = 0;
  for (std::map<double,double>::reverse_iterator i = distribution.rbegin(); i != distribution.rend(); ++i, aux_counter++) {
    values[aux_counter] = ((*i).second);
  }

  aux_counter = 0;
  for (std::map<double,double>::iterator i = distribution.begin(); i != distribution.end(); ++i, aux_counter++) {
    keys[aux_counter] = ((*i).first);
  }

  std::partial_sum(values.begin(), values.end(), std::back_inserter(counts_partial_sums));


  for (size_t i = 0; i < thresholds.size(); i++) {
    if (thresholds[i] < cutoff) {
      p_values[i] = 1.0;
      continue;
    }
    //NOTE: it seems that cache doesn't make any speedup at all: the same speed with it or without on the same tests so I have deleted it.
    // It's reasonable: complexity of search in big std::map is higher than followin operations. 
    double counts_sum = 0.0;

    std::vector<double>::iterator lower = std::lower_bound(keys.begin(), keys.end(), thresholds[i] * DISCRETIZATION_VALUE);
              
    counts_sum = counts_partial_sums[std::distance(lower, keys.end()) - 1];
    
    p_values[i] = counts_sum/volume;

  }
  return;
}





