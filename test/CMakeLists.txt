# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

project (edda_test)

include_directories(${CMAKE_SOURCE_DIR}/src)

set(TEST_SOURCES
    test_gaussian.cpp
    test_variant.cpp
    test_thrust.cpp
    )

foreach (test ${TEST_SOURCES})
  get_filename_component(tname ${test} NAME_WE)
  add_executable(${tname} ${test})
  target_link_libraries(${tname} ${LINK_LIBS})

  add_test(NAME ${tname}
    COMMAND ${tname}
    )
endforeach (test)


if (EDDA_WITH_VTK)
  set(sources
    gmm_raw2vtk
    test_vtk_pipeline
    #test_vtkEddaReader
  )

  IF(APPLE)
    # The primitive VTK renderer cannot catch keyboard events
    # without packed into the mac bundle
    SET(EXECUTABLE_FLAG MACOSX_BUNDLE)
  ENDIF(APPLE)

  foreach (file ${sources})
    get_filename_component(name ${file} NAME_WE)
    add_executable ( ${name} ${EXECUTABLE_FLAG} ${file} )
    target_link_libraries ( ${name} ${LINK_LIBS} ${VTK_LIBRARIES})
  endforeach ()

endif ()

