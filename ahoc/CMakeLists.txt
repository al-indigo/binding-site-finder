set(ahocorasick_C_SRC ahocorasick.c node.c)
add_library(ahoc_C_LIB STATIC ${ahocorasick_C_SRC})

set(ahocorasick_CPP_SRC AhoCorasickPlus.cpp)
add_library(ahoc_CPP_LIB STATIC ${ahocorasick_CPP_SRC})
target_link_libraries(ahoc_CPP_LIB ahoc_C_LIB)