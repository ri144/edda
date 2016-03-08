# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

set(sources
  uncertainIsocontour.cu
  #streamlineGaussian.cu
)

include_directories(${CMAKE_SOURCE_DIR}/src)

# Force compile .cu files
add_definitions("-x c++")  


foreach (file ${sources})

  get_filename_component(name ${file} NAME_WE)

  if (EDDA_BUILD_WITH_CUDA)
    cuda_add_executable ( ${name} ${file} )
  else ()
    SET_SOURCE_FILES_PROPERTIES(${file} PROPERTIES LANGUAGE CXX)   
    add_executable ( ${name} ${file} )
  endif ()

  target_link_libraries ( ${name} edda ${LINK_LIBS} )
endforeach ()

if (EDDA_WITH_VTK)
  add_subdirectory (vtk)
endif ()

