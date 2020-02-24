cmake_minimum_required(VERSION 3.6)
project("Neighborhood search")

add_library(neighborhood_search src/neighborhood_search.c)
target_include_directories(neighborhood_search PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(neighborhood_search
                      PUBLIC bov pthread)

add_library(neighborhood_search_for_mac src/neighborhood_search_for_mac.c)
target_include_directories(neighborhood_search_for_mac PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(neighborhood_search_for_mac
                      PUBLIC bov)