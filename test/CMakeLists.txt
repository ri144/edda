# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

project (edda_test)

include_directories(${CMAKE_SOURCE_DIR}/src)

# Not CUDA
if (NOT EDDA_BUILD_WITH_CUDA)
  set(TEST_SOURCES
      test_gaussian.cpp
      test_thrust.cpp
      test_NdArray.cpp
  )

  foreach (test ${TEST_SOURCES})
    get_filename_component(tname ${test} NAME_WE)
    add_executable(${tname} ${test})
    target_link_libraries(${tname} edda ${LINK_LIBS})

    add_test(NAME ${tname}
      COMMAND ${tname}
      )
  endforeach (test)
endif ()

# CUDA
if (EDDA_BUILD_WITH_CUDA)
  set(TEST_SOURCES_CUDA
      test_gaussian.cu
      test_thrust.cu
      )

  foreach (test ${TEST_SOURCES_CUDA})
    get_filename_component(tname ${test} NAME_WE)
    cuda_add_executable(${tname} ${test})
    target_link_libraries(${tname} edda ${LINK_LIBS})

    add_test(NAME ${tname}
      COMMAND ${tname}
      )
  endforeach (test)
endif ()

# VTK
if (EDDA_WITH_VTK)
  set(sources
    gmm_raw2vtk.cpp
    #test_vtk_pipeline.cpp
  )

  IF(APPLE)
    # The primitive VTK renderer cannot catch keyboard events
    # without packed into the mac bundle
    SET(EXECUTABLE_FLAG MACOSX_BUNDLE)
  ENDIF(APPLE)

  foreach (file ${sources})
    get_filename_component(name ${file} NAME_WE)

    if (EDDA_BUILD_WITH_CUDA)
      cuda_add_executable ( ${name} ${EXECUTABLE_FLAG} ${file} )
    else ()
      add_executable ( ${name} ${EXECUTABLE_FLAG} ${file} )
    endif ()

    target_link_libraries ( ${name} edda ${LINK_LIBS} )
  endforeach ()

endif ()

