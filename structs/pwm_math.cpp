#include "pwm_math.h"

#include <iostream>
#include <cmath>
#include <numeric>
#include <utility>
#include <stdlib.h>
#include <algorithm>

/* TODO: NOTE: I know that using doubles in comparisons (or find) is not semantically right, need to fix sometime */

//Approximated inverse Gauss error function
inline double inverf (double x) {
  //NOTE: don't want to use c++11 stdrt just for this.
  int sign = (x > 0) ? 1 : ((x < 0) ? -1 : 0);
  x = (x > 0) ? x : -x;
  double a = 0.147; //from wiki: this is more accurate value for approximation
  double tmp = (2.0 / (M_PI * a) + (log(1.0 - x * x)) / 2.0);
  double part0 = tmp * tmp;
  double part = -2.0 / (M_PI * a) - log(1.0 - x * x) / 2.0 + sqrt(-1.0 / a * log(1.0 - x * x) + part0);
//  std::cout << x << " " << -2.0 / (M_PI * a) << " " << - log(1.0 - x * x) / 2.0 << " " << sqrt(-1.0 / a * log(1.0 - x * x) + part0) << std::endl;
//  std::cout << "sign: " << sign << " a: " << a << " tmp: " << tmp << " part0: " << part0 << " part: " << part << " PI: " << M_PI << std::endl;

  return sign * sqrt(part);
}

inline static std::map<double, double> recalc_score_hash (std::map<double, double> scores, pwmMatrix& matrix, size_t row, double least_sufficient) {
  std::map<double, double> new_scores;
  
  for (std::map <double, double>::iterator i = scores.begin(); i != scores.end(); ++i) {
    double score = (*i).first;
    double count = (*i).second;
    for (int k = 0; k < matrix.cols; k++) {
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

inline static std::map<double, double> count_distribution_after_threshold (pwmMatrix& matrix, double * scores_optimistic, double threshold) {
  std::map<double, double> scores;
  scores.insert(std::make_pair(0.0, 1.0));
  for (int row = 0; row < matrix.rows; row++) {
    scores = recalc_score_hash(scores, matrix, row, threshold - scores_optimistic[row]);
  }
  return scores;
}

inline static double score_mean (pwmMatrix& matrix) {
  double result = 0.0;
  for (int row = 0; row < matrix.rows; row++) {
    double row_sum = 0.0;
    for ( int col = 0; col < matrix.cols; col++) {
      row_sum += matrix(row,col);
    }
    result += row_sum/4.0;
  }
  return result;
}

inline static double score_variance (pwmMatrix& matrix) {
  double result = 0.0;
  for (int row = 0; row < matrix.rows; row++) {
    double mean_value = 0.0;
    double mean_square = 0.0;
    
    for ( int col = 0; col < matrix.cols; col++) {
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
//  std::cout << "sigma: " << sigma <<std::endl;
//  std::cout << "test inverf param: " << 1.0 - 2.0 * pvalue << std::endl;
  double n_ = inverf(1.0 - 2.0 * pvalue)*sqrt(2.0);
//  std::cout << "inverf: " << n_ << std::endl;
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
      cnt_distribution = count_distribution_after_threshold(matrix, scores_optimistic, approximate_threshold);
      max_pvalue *= 2; // if estimation counted too small amount of words - try to lower threshold estimation by doubling pvalue
      
      values_sum = 0.0;
      for (std::map<double, double>::iterator i = cnt_distribution.begin(); i != cnt_distribution.end(); ++i) {
        values_sum += (*i).second;
      }
//      std::cout << "values sum: " << values_sum << ", expected:" << look_for_count << ", approx threshold=" << approximate_threshold << std::endl;
    }

    return cnt_distribution;
}

double threshold_by_pvalue (double p_value, double * scores_optimistic, pwmMatrix& matrix) {
  double threshold = 0.0;
  std::map<double, double> scores_hash = count_distribution_under_pvalue(p_value, scores_optimistic, matrix);
  size_t total_count = 0;
  std::vector<double> values, keys, scores_partial_sums;
//  values.reserve(scores_hash.size());
  scores_partial_sums.reserve(scores_hash.size());
  for (std::map<double,double>::iterator i = scores_hash.end(); i != scores_hash.begin(); --i) {
//    std::cout << (*i).first << "   " << (*i).second << std::endl;
    total_count += (*i).second;
    values.push_back((*i).second);
    keys.push_back((*i).first);
  }
  std::vector<double>::iterator test = std::partial_sum<std::vector<double>::iterator, std::vector<double>::iterator>(values.begin(), values.end(), scores_partial_sums.end());
  std::cout << "Size of hashmap is: " << scores_hash.size() << std::endl;
/*  for (test; test != scores_partial_sums.begin(); --test) {
    std::cout << *test << " " ;
  }*/
  for (std::vector<double>::iterator i = scores_partial_sums.begin(); i != test; i++) {
//    std::cout << (double)*i << " " << i - scores_partial_sums.begin() << std::endl;
  }
  
  std::pair<std::vector<double>::iterator, std::vector<double>::iterator> bounds;
  bounds = std::equal_range(scores_partial_sums.begin(), test, floor(p_value * vocabularyVolume(matrix)));
  
//   std::cout << "Need P-Value: " << p_value << " with " << floor(p_value*vocabularyVolume(matrix)) << " words" << std::endl;
//   std::cout << "lower bound: " << (bounds.first) - scores_partial_sums.begin() << ", upper bound: " << (bounds.second) - scores_partial_sums.begin() << std::endl;
//   std::cout << "lower distance: " << (size_t)abs(std::distance(bounds.first, scores_partial_sums.begin())) << ", upper distance: " << (size_t)abs(std::distance(bounds.second, scores_partial_sums.begin())) << std::endl;
  
  size_t l_d = (size_t)abs(std::distance(bounds.first, scores_partial_sums.begin()));
  size_t u_d = (size_t)abs(std::distance(bounds.second, scores_partial_sums.begin()));
  
//   std::cout << "lower bound: " << (*bounds.first) << ", upper bound: " << (*bounds.second) << std::endl;
  
  std::map<double,double>::iterator it = scores_hash.begin();
//   std::cout << "      threshold: "<< " " << (keys[l_d])/DISCRETIZATION_VALUE << std::endl;
//   std::cout << "      actual p-value is: " << (double)scores_partial_sums[l_d]/(double)vocabularyVolume(matrix) << std::endl;

  
  threshold = (keys[l_d])/DISCRETIZATION_VALUE;
  
  return threshold;
}

std::vector<double> pvalues_by_thresholds(std::vector<double>& thresholds, double * scores_optimistic, pwmMatrix& matrix) {
  double min_threshold = *std::min_element(thresholds.begin(), thresholds.end());
  std::map<double, double> counts = count_distribution_after_threshold(matrix, scores_optimistic, min_threshold*DISCRETIZATION_VALUE);
  
  double volume = (double)vocabularyVolume(matrix);
  std::vector<double> p_values;
  p_values.reserve(thresholds.size());
  
  std::map<double, double> cache;
  
  for (std::vector<double>::iterator i = thresholds.begin(); i != thresholds.end(); ++i) {
    std::map<double, double>::iterator search = cache.find(*i);
    if (search != cache.end()) {
      p_values.push_back(search->second);
    } else {
      double counts_sum = 0.0;
      std::map<double, double>::iterator lower = counts.lower_bound(*i * DISCRETIZATION_VALUE);
      for (lower; lower != counts.end(); lower++) {
        counts_sum += (*lower).second;
      }
      p_values.push_back(counts_sum/volume);
      cache.insert(std::make_pair(*i, counts_sum/volume));
    }
  }
  return p_values;
}






