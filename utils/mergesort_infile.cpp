#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

void merge_sort(std::string filename1in, std::string filename2in, std::string filename_out) {
    // Open Files
    std::ifstream fin(filename1in);
    std::ifstream fin2(filename2in);
    std::ofstream fout(filename_out);

    std::string line;                // string to hold line from file
//    size_t i = 1;                  // line counter
    size_t in1 = 0;
    size_t in2 = 0;
    if(fin) {
        getline(fin,line);
        std::stringstream ss(line);
        ss >> in1;
    }
    if(fin2) {
        getline(fin2,line);
        std::stringstream ss(line);
        ss >> in2;
    }
    bool first = true;          // bool to catch when a file closes
    while(fin || fin2) {
        if(fin && fin2) {
            if(in2 <= in1) {
                fout << in2 << std::endl;
                getline(fin2, line);
                std::stringstream ss(line);
                ss >> in2;
            }
            else {
                fout << in1 << std::endl;
                getline(fin, line);
                std::stringstream ss(line);
                ss >> in1;
            }
        } else {
            if(first) {     // first time through the else
                if(!fin) {
                  fout << in2 << std::endl;
                } else if(!fin2) {
                  fout << in1 << std::endl;
                }
            } else {
                fout << line << std::endl;
            }

            // get the next line from the file that is open
            if(fin) {
              getline(fin, line);
            } else if(fin2) {
              getline(fin2, line);
            }
            first = false;  // only print line from now on, don't print in1 or in2
        }
    }
}
