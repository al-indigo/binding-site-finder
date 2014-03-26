cmake_minimum_required(VERSION 2.8)

add_executable(benchmark bench.cpp)

add_custom_target(copy)
add_custom_command(TARGET copy PRE_BUILD                              # Adds a post-build event to MyTest
                   COMMAND ${CMAKE_COMMAND} -E copy         # which executes "cmake - E copy_if_different..."
                   "${PROJECT_SOURCE_DIR}/benchmark/benchmark.json"      # <--this is in-file
                    ${PROJECT_BINARY_DIR}/benchmark/benchmark.json)         # <--this is out-file path

add_dependencies(benchmark copy ahoc_CPP_LIB structs_LIB utils_CPP_LIB)



target_link_libraries(benchmark ahoc_CPP_LIB structs_LIB utils_CPP_LIB functions_LIB)
