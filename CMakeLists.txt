# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

cmake_minimum_required (VERSION 2.8)

# Set the project's name.
project (EDDA)

# Tell the compiler to use C++11.
if (UNIX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
endif ()

# Set libraries' output path to lib.
set (LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib")

# Set variables for platform detection.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set (OS_LINUX ON)
endif ()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set (OS_MACOSX ON)
endif ()

if (WIN32)
  set (OS_WIN ON)
endif ()

# Boost
find_package(Boost REQUIRED)

# Provide options that the user can optionally select.
option (BUILD_SHARED_LIBS "Build into shared libraries." ON)

option (EDDA_BUILD_TESTS "Build the test programs." OFF)
option (EDDA_BUILD_EXAMPLES "Build the example programs." OFF)

# Make the configuration file.
configure_file ("${PROJECT_SOURCE_DIR}/edda.h.in"
                "${PROJECT_SOURCE_DIR}/src/edda.h")

# Add include directories to the build.
include_directories ("${PROJECT_SOURCE_DIR}/src")

# Add the src subdirectory to the build.
add_subdirectory (src)

# If EDDA_BUILD_TESTS is on, enable testing.
if (EDDA_BUILD_TESTS)
  enable_testing ()
  add_subdirectory (test)
endif ()

# If EDDA_BUILD_EXAMPLES is on, add the examples subdirectory to the build.
if (EDDA_BUILD_EXAMPLES)
  add_subdirectory (examples)
endif ()
