# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

project (edda_test)

include_directories(${CMAKE_SOURCE_DIR}/src)
add_definitions(-DSAMPLE_DATA_PATH="${CMAKE_SOURCE_DIR}/sample_data")


if (NOT EDDA_BUILD_WITH_CUDA)
  set(TEST_SOURCES
    test_gaussian.cpp
    test_thrust.cpp
    test_NdArray.cpp
  )
else ()
  # CUDA
  set(TEST_SOURCES
      test_gaussian_cuda.cu
      test_thrust_cuda.cu
      test_NdArray_cuda.cu
      )
endif ()

# VTK
if (EDDA_WITH_VTK)
  set(TEST_SOURCES ${TEST_SOURCES}
      test_vtk_filters.cpp
  )
endif()

####
foreach (test ${TEST_SOURCES})
  get_filename_component(tname ${test} NAME_WE)
  if (NOT EDDA_BUILD_WITH_CUDA)
    add_executable(${tname} ${test})
  else ()
    cuda_add_executable(${tname} ${test})
  endif ()
  target_link_libraries(${tname} edda ${LINK_LIBS})

  add_test(NAME ${tname}
    COMMAND ${tname}
    )
endforeach (test)

