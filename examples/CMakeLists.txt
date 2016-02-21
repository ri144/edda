# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

set(sources
  uncertainIsocontour
  streamlineGaussian
)

include_directories(${CMAKE_SOURCE_DIR}/src)

foreach (file ${sources})
  get_filename_component(name ${file} NAME_WE)
  add_executable ( ${name} ${file} )
  target_link_libraries ( ${name} ${LINK_LIBS} )
endforeach ()

if (EDDA_WITH_VTK)
  add_subdirectory (vtk)
endif ()
