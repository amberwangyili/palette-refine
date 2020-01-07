cmake_minimum_required(VERSION 3.12)
project(palette_refine)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_POSTFIX "${CMAKE_RELEASE_POSTFIX}")

find_package(OpenCV)

include_directories(${OpenCV_INCLUDE_DIRS})

# Judging platform environment
if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "Windows settings are processed here!")
else()
    message(STATUS "Linux/Unix compiler settings are processed here!")
endif()

#Output Messages for debug the Cmake
message(STATUS "operation system is : ${CMAKE_SYSTEM}")
message(STATUS "current platform is : ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMake version is : ${CMAKE_SYSTEM_VERSION}")
message(STATUS "The program main directory is : ${ROOT}")

# Current compiler
message(STATUS "C compiler is : ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler is : ${CMAKE_CXX_COMPILER}")



add_executable(palette_refine main.cpp vec3.h nearestPoint.h cxxopt.h)

target_link_libraries(palette_refine ${OpenCV_LIBS})

target_link_libraries(palette_refine -lnlopt -lm -lstdc++)