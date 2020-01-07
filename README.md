```cmake
# Minimum requirement of CMake version : 3.0.0
cmake_minimum_required(VERSION 3.0.0)

#Project name and version number
project(Project VERSION 0.1.0)

# Set the ROOT and subdirectory
# you should put the CMakeList.txt in these file directories
set(ROOT ${PROJECT_SOURCE_DIR})
# 可执行文件输出路径
set(EXECUTABLE_OUTPUT_PATH ${ROOT}/Bin)
# 可执行文件输出路径，和上语句作用类似
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT}/Bin)
# 生成的库文件输出路径
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT}/Lib)
# 主要针对VS或者XCode之类的IDE,指定库文件路径
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT}/Lib)
# 如果是debug模式，加后缀_d
set(CMAKE_DEBUG_POSTFIX "_d")
# 如果是release模式，加后缀_r
set(CMAKE_RELEASE_POSTFIX "_r")

#Choose different compilation configurations according to VS compilation
set(CMAKE_BUILD_POSTFIX "${CMAKE_RELEASE_POSTFIX}") 

# For lib source files
add_subdirectory(${ROOT}/Src/Libs ${ROOT}/build/Libs)
# For main file of the function apps
add_subdirectory(${ROOT}/Src/Apps ${ROOT}/build/Apps)

# Judging platform environment
if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "Windows settings are processed here!")
else()
    message(STATUS "Linux/Unix compiler settings are processed here!") 
    # Turn on c++ 11
    set(CMAKE_CXX_STANDARD 11)  
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    # Compile default settings, more settings can be selected in the cmake command
    add_compile_options(-w)
    add_compile_options(-m64)
    add_compile_options(-lz)
    add_compile_options(-std=c++11)
endif()

#Output Messages for debug the Cmake
message(STATUS "operation system is : ${CMAKE_SYSTEM}") 
message(STATUS "current platform is : ${CMAKE_SYSTEM_NAME}")        
message(STATUS "CMake version is : ${CMAKE_SYSTEM_VERSION}") 
message(STATUS "The program main directory is : ${ROOT}")

# Current compiler
message(STATUS "C compiler is : ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler is : ${CMAKE_CXX_COMPILER}")   
```