cmake_minimum_required(VERSION 3.29)
project(Quadratic_Zeeman_Effect C)

set(CMAKE_C_STANDARD 11)

include_directories(Code)

add_executable(Quadratic_Zeeman_Effect
    Code/main.c
    Code/perturbation_equation.c
    Code/perturbation_equation.h)
