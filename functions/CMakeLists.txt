set(functions_SRC init_ahoc.cpp ahoc_search.cpp merge_files.cpp recalc_scores_p_values.cpp format_output.cpp write_status.cpp)
add_library(functions_LIB STATIC ${functions_SRC})
target_link_libraries(functions_LIB ahoc_CPP_LIB structs_LIB utils_CPP_LIB)
