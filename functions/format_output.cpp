#include "functions.h"

#include <iostream>
#include <ostream>
#include <cmath>
#include <vector>
#include <cstdio>


void format_bed(std::ostream& fout,
                std::string& chromoname,
                std::vector<size_t>& positions_to_read_again,
                std::vector<double>& pvaluesFw,
                std::vector<double>& pvaluesRev,
                std::vector<double>& scoresFw,
                std::vector<double>& scoresRev) {
  std::vector<size_t>::iterator positerator = positions_to_read_again.begin()++;

  double tstart, tstop, ttime;
  tstart = (double)clock()/CLOCKS_PER_SEC;
  for (int k = 0; k < positions_to_read_again.size() ; k++) {
    double pv = std::min(pvaluesFw[k], pvaluesRev[k]);

    int narrowPeakScore = 0;
    if (pv > 0.001) {
      narrowPeakScore = 100; 
    } else if (pv < 0.00001) {
      narrowPeakScore = 1000;
    } else {
      narrowPeakScore = (int) (100.0 + 900.0 * (-3.0 - log10f(pv)) / 2.0);
    }
    
    fout      << chromoname << "\t" 
              << *positerator << "\t" << *positerator << "\t"
              << ".\t"
              << narrowPeakScore << "\t" 
              << (scoresFw[k] > scoresRev[k] ? "+" : "-") << "\t" 
              << std::max(scoresFw[k], scoresRev[k]) << "\t" 
              << pv << "\t" << "-1\t-1" << "\n";
    positerator++;
  }
  tstop = (double)clock()/CLOCKS_PER_SEC;
  std::cout << "Formatted out for " << tstop - tstart << " seconds" << std::endl; 
}

void format_bed(FILE * fout,
                std::string& chromoname,
                std::vector<size_t>& positions_to_read_again,
                std::vector<double>& pvaluesFw,
                std::vector<double>& pvaluesRev,
                std::vector<double>& scoresFw,
                std::vector<double>& scoresRev) {
  std::vector<size_t>::iterator positerator = positions_to_read_again.begin()++;

  double tstart, tstop, ttime;
  tstart = (double)clock()/CLOCKS_PER_SEC;
  for (int k = 0; k < positions_to_read_again.size() ; k++) {
    double pv = std::min(pvaluesFw[k], pvaluesRev[k]);

    int narrowPeakScore = 0;
    if (pv > 0.001) {
      narrowPeakScore = 100; 
    } else if (pv < 0.00001) {
      narrowPeakScore = 1000;
    } else {
      narrowPeakScore = (int) (100.0 + 900.0 * (-3.0 - log10f(pv)) / 2.0);
    }
    fprintf(fout, "%s\t%lu\t%lu\t.\t%u\t%s\t%.6f\t%.8f\t-1\t-1\n", chromoname.c_str(), *positerator, *positerator, narrowPeakScore, (scoresFw[k] > scoresRev[k] ? "+" : "-"), std::max(scoresFw[k], scoresRev[k]), pv);
              
    positerator++;
  }
  tstop = (double)clock()/CLOCKS_PER_SEC;
  std::cout << "Formatted out for " << tstop - tstart << " seconds" << std::endl; 
}