cmake_minimum_required(VERSION 2.8)

add_executable(benchmark bench.cpp)

add_custom_target(copy)
add_custom_command(TARGET copy PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                   "${PROJECT_SOURCE_DIR}/benchmark/benchmark.json"
                    ${PROJECT_BINARY_DIR}/benchmark/benchmark.json)

add_dependencies(benchmark copy ahoc_CPP_LIB structs_LIB utils_CPP_LIB)


IF(APPLE)
target_link_libraries(benchmark ahoc_CPP_LIB structs_LIB utils_CPP_LIB functions_LIB)
ELSE()
target_link_libraries(benchmark ahoc_CPP_LIB structs_LIB utils_CPP_LIB functions_LIB rt)
ENDIF(APPLE)