#include "mergesort_infile.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cstdio>
#include <iomanip>

void merge_sort(std::string filename1in, std::string filename2in, std::string filename_out) {

    std::ifstream fin1(filename1in.c_str());
    std::ifstream fin2(filename2in.c_str());
    std::ofstream fout(filename_out.c_str());
    //NOTE: working with HEX because it takes much less space and IO is critical
    fin1 >> std::setbase(16);
    fin2 >> std::setbase(16);
    fout << std::setbase(16);

    std::merge(std::istream_iterator<size_t>(fin1),
               std::istream_iterator<size_t>(),
               std::istream_iterator<size_t>(fin2),
               std::istream_iterator<size_t>(),
               std::ostream_iterator<size_t>(fout, "\n"));
    fin1.close();
    fin2.close();
    fout.close();
    remove(filename1in.c_str());
    remove(filename2in.c_str());
}