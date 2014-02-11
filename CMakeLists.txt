cmake_minimum_required(VERSION 2.8)
project(binding-site-finder)
set (PROJECT binding-site-finder)

#set (Boost_USE_STATIC_LIBS ON)
#set (Boost_USE_MULTITHREADED ON)
#find_package (Boost COMPONENTS regex REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(ahoc)
add_subdirectory(structs)
add_subdirectory(utils)

add_executable(${PROJECT} main.cpp)
add_dependencies(${PROJECT} ahoc_CPP_LIB structs_LIB utils_CPP_LIB)

target_link_libraries(${PROJECT} ahoc_CPP_LIB structs_LIB utils_CPP_LIB)

#target_link_libraries(${PROJECT} ${Boost_LIBRARIES})
