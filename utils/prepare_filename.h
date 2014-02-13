#ifndef PREPARE_FILENAME_H
#define PREPARE_FILENAME_H

#include <string>
#include <time.h>


std::string prepare_filename (std::string prefix, std::string suffix, int number, time_t * time = NULL);

#endif