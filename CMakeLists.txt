cmake_minimum_required(VERSION 3.22)
project(kp_da)

set(CMAKE_CXX_STANDARD 20)
set(GCC_FLAGS "-O2 -pedantic -std=c++20 -Wall -Wextra -lm")
set(CMAKE_CXX_FLAGS "${GCC_FLAGS}")

add_executable(kp_da main.cpp utils.h astar.h)
