cmake_minimum_required(VERSION 3.14)

project(single-matrix-affine VERSION 1.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)

#
# SETTINGS
#

# List of valid build types
set(VALID_BUILD_TYPES Debug Release RelWithDebInfo MinSizeRel)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as no build type was specified")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Add warnings for debug builds
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wuninitialized -Wunused -Wshadow -Wconversion")
    message(STATUS "Building in debug mode")
endif()

# Check if the provided CMAKE_BUILD_TYPE is in the valid list
if(NOT CMAKE_BUILD_TYPE IN_LIST VALID_BUILD_TYPES)
    message(FATAL_ERROR "Invalid build type: ${CMAKE_BUILD_TYPE}. Valid options are: ${VALID_BUILD_TYPES}")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Add source files and executable.
file(GLOB SOURCES src/*.cpp src/*/*.cpp)
add_executable(benchmark ${SOURCES})

