set(utils_CPP_SRC mergesort_infile.cpp)
add_library(utils_CPP_LIB STATIC ${utils_CPP_SRC})
target_link_libraries(utils_CPP_LIB)