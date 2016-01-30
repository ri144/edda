# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

set(SRC_FILES "")
list ( APPEND SRC_FILES
    vtkEddaReader.cpp
    vtkSamplingArray.cpp
    )

file (GLOB_RECURSE HEADER_FILES "*.h")

set(PROJECT_NAME_VTK eddaVTK)

add_library (${PROJECT_NAME_VTK} ${LIB_TYPE} ${SRC_FILES} ${HEADER_FILES})

set_target_properties (${PROJECT_NAME_VTK} PROPERTIES
    DEBUG_OUTPUT_NAME "${PROJECT_NAME_VTK}_d"
    RELEASE_OUTPUT_NAME "${PROJECT_NAME_VTK}")

IF(APPLE)
  # The primitive VTK renderer cannot catch keyboard events
  # without packed into the mac bundle
  SET(EXECUTABLE_FLAG MACOSX_BUNDLE)
ENDIF(APPLE)

#----------------------------------------
# test
#----------------------------------------
#add_executable(test_vtkEddaReader test_vtkEddaReader.cpp )
#target_link_libraries( test_vtkEddaReader ${PROJECT_NAME} ${PROJECT_NAME_VTK} ${VTK_LIBRARIES} )
add_executable(test_vtkEddaReader ${EXECUTABLE_FLAG} test_vtkEddaReader.cpp ${SRC_FILES})
target_link_libraries( test_vtkEddaReader ${PROJECT_NAME} ${VTK_LIBRARIES} )

add_executable(test_pipeline ${EXECUTABLE_FLAG} test_pipeline.cpp ${SRC_FILES})
target_link_libraries( test_pipeline ${PROJECT_NAME} ${VTK_LIBRARIES} )

install (
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION include/edda/vtk
    FILES_MATCHING PATTERN "*.h")

#install (
#    TARGETS ${PROJECT_NAME_VTK}
#    DESTINATION lib/edda)


if (${EDDA_BUILD_PARAVIEW_PLUGINS})
  ADD_PARAVIEW_PLUGIN(eddaReader "1.0"
    SERVER_MANAGER_XML vtkEddaReader.xml
    SERVER_MANAGER_SOURCES ${SRC_FILES})
  TARGET_LINK_LIBRARIES(eddaReader LINK_PRIVATE ${PROJECT_NAME})
endif ()

