cmake_minimum_required(VERSION 3.24)
project(math_lab_1 C)

set(CMAKE_C_STANDARD 11)

add_executable(math_lab_1 main.c spline.c cmath.h cmathmsg.c quanc8.c)
target_link_libraries(math_lab_1 m)