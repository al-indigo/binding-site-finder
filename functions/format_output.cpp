#include "functions.h"

#include <iostream>
#include <ostream>
#include <cmath>
#include <vector>

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
  //          narrowPeakScore = 100 + (int) (900 * ((0.001 - pv)/(0.00099) ) );
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