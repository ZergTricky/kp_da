cmake_minimum_required(VERSION 3.22)
project(kp_da)

set(CMAKE_CXX_STANDARD 20)
set(GCC_FLAGS "-Os -pedantic -pg -std=c++20 -Wall -Wextra -lm")
set(CMAKE_CXX_FLAGS "${GCC_FLAGS}")

add_executable(solution main.cpp search.h utils.h)
add_executable(gen gen.cpp)
add_executable(validator validator.cpp)