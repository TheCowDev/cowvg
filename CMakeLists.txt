cmake_minimum_required(VERSION 3.25)
project(cowvg C)

set(CMAKE_C_STANDARD 11)

include_directories(.)
include_directories(cowvg)

# Set C flags for Release mode
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")

add_executable(cowvg
        cowvg/main.c
        cowvg/cowvg.h cowvg/cowvg.c)
