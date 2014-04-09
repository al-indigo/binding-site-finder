cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS                "-Wall -std=c++11 -msse3 -msse4 -msse4.1 -mavx -march=native")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

project(bsf)
set (PROJECT bsf)

add_subdirectory(ahoc)
add_subdirectory(structs)
add_subdirectory(utils)
add_subdirectory(benchmark)
add_subdirectory(functions)

add_executable(${PROJECT} bsf.cpp)
add_dependencies(${PROJECT} ahoc_CPP_LIB structs_LIB utils_CPP_LIB functions_LIB)

target_link_libraries(${PROJECT} ahoc_CPP_LIB structs_LIB utils_CPP_LIB functions_LIB)