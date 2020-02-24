cmake_minimum_required(VERSION 3.6)
project("Neighborhood search")

add_library(neighborhood_search src/neighborhood_search.c)
target_include_directories(neighborhood_search PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(neighborhood_search PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(neighborhood_search
                      PUBLIC bov pthread)