cmake_minimum_required(VERSION 2.8)

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