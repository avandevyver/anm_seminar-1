cmake_minimum_required(VERSION 3.0)
project(Pthread)

add_library(pthread src/pthread.c)
target_include_directories(pthread PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(pthread PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)