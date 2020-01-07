cmake_minimum_required(VERSION 3.12)
project(palette_refine)

set(CMAKE_CXX_STANDARD 14)
#find_library(OpenCV)
find_package(OpenCV)

include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(palette_refine main.cpp vec3.h nearestPoint.h cxxopt.h)

target_link_libraries(palette_refine ${OpenCV_LIBS})

target_link_libraries(palette_refine -lnlopt -lm -lstdc++)