cmake_minimum_required(VERSION 2.8)

project(bsf)
set (PROJECT bsf)

if(MSVC AND (MSVC10 OR MSVC11 OR MSVC12))
  message(FATAL_ERROR "Sorry, this program can not be built using Visual Studio (cygwin/mingw are also untested)")
endif()



include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
  set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
  set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")


  add_subdirectory(ahoc)
  add_subdirectory(structs)
  add_subdirectory(utils)
  #add_subdirectory(benchmark)
  add_subdirectory(functions)

  add_executable(${PROJECT} bsf.cpp)
  add_dependencies(${PROJECT} ahoc_CPP_LIB structs_LIB utils_CPP_LIB functions_LIB)

  target_link_libraries(${PROJECT} ahoc_CPP_LIB structs_LIB utils_CPP_LIB functions_LIB pthread)
  
else()
  message(FATAL_ERROR "You should have a compiler with C++11 standard support to compile bsf\n(gcc-4.7+, clang-3.2+, Intel C compiler 13.0+)\nNow you are trying to use ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
  
endif()