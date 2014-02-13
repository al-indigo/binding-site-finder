#include "prepare_filename.h"
#include <sstream>

std::string prepare_filename (std::string prefix, std::string suffix, int number, time_t * t) {
  std::stringstream temp_ss_stream;
  if (t == NULL) {
    temp_ss_stream << prefix << suffix << time(NULL) << "-" << number;
  } else {
    temp_ss_stream << prefix << suffix << *t << "-" << number;
  }
  std::string tempfilename;
  temp_ss_stream >> tempfilename;
  return tempfilename;
}