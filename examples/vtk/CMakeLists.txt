# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

set(sources
  uncertainIsocontourVTK.cpp
  streamlineMonteCarlo.cpp
  streamlineGaussianVTK.cpp
)

IF(APPLE)
  # The primitive VTK renderer cannot catch keyboard events
  # without packed into the mac bundle
  SET(EXECUTABLE_FLAG MACOSX_BUNDLE)
ENDIF(APPLE)

foreach (file ${sources})
  get_filename_component(name ${file} NAME_WE)
  add_executable ( ${name} ${EXECUTABLE_FLAG} ${file} )
  target_link_libraries ( ${name} ${PROJECT_NAME} ${VTK_LIBRARIES})
endforeach ()

