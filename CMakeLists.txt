# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

cmake_minimum_required (VERSION 2.8)

# Set the project's name.
project (edda)

### Provide options that the user can optionally select.
option (EDDA_BUILD_TESTS "Build the test programs." OFF)
option (EDDA_BUILD_EXAMPLES "Build the example programs." ON)
option (EDDA_BUILD_SHARED_LIBS "Build into shared libraries." ON)
option (EDDA_BUILD_PARAVIEW_PLUGINS "Build Paraview Plugins." OFF)
option (EDDA_WITH_VTK "Build with VTK integration." OFF)
option (EDDA_WITH_OPENMP "Build with OpenMP parallelization." ON)
###

list(APPEND LINK_LIBS ${PROJECT_NAME})

# Tell the compiler to use C++11.
if (UNIX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall)
endif ()

if (EDDA_WITH_OPENMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif ()

## Thrust
set(THRUST_INCLUDE_DIR "/usr/local/include" CACHE PATH "Thrust root dir")
INCLUDE_DIRECTORIES(${THRUST_INCLUDE_DIR})

## Thrust with TBB backend
option (THRUST_DEVICE_SYSTEM_TBB "Build with Cuda Thrust parallelization." OFF)
if (THRUST_DEVICE_SYSTEM_TBB)
  set(TBB_INSTALL_PATH "/usr/local" CACHE PATH "TBB root dir")
  INCLUDE_DIRECTORIES(${TBB_INSTALL_PATH}/include)
  LINK_DIRECTORIES(${TBB_INSTALL_PATH}/lib)
  add_definitions(-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB)
  list(APPEND LINK_LIBS tbb)
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

if (${EDDA_BUILD_SHARED_LIBS})
    set(LIB_TYPE SHARED)
else ()
    set(LIB_TYPE STATIC)
endif ()

# Load settings for the boost project.
find_package (Boost REQUIRED)
if (Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
  link_libraries (${Boost_LIBRARIES})
endif ()


# Choose either Paraview or VTK library
if (EDDA_BUILD_PARAVIEW_PLUGINS)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})

else ()
  # VTK integration
  if (EDDA_WITH_VTK )
    find_package(VTK)
    include (${VTK_USE_FILE})
  endif ()
endif ()

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

