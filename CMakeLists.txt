cmake_minimum_required(VERSION 3.0.0)
project(glowstone VERSION 0.1.0)

add_library(glowstone source/glowstone.cpp source/more_math.cpp)
target_include_directories(glowstone PUBLIC include)

add_executable(main_test test/main_test.cpp)
target_link_libraries(main_test glowstone)
