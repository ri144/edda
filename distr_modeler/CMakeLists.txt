# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

project (edda_distr_modeler)

include_directories(${CMAKE_SOURCE_DIR}/src)
add_definitions(-DSAMPLE_DATA_PATH="${CMAKE_SOURCE_DIR}/sample_data")

##add VTK to the project
#find_package(VTK)
#include (${VTK_USE_FILE})
#message("vtk location: " ${VTK_USE_FILE})

find_package(NetCDF)

set(CONVERTER_SOURCES ensembleExample.cpp spatialDecompositionExample.cpp)

####
foreach (ldsrc ${CONVERTER_SOURCES})
  get_filename_component(fname ${ldsrc} NAME_WE)
  if (NOT EDDA_BUILD_WITH_CUDA)
    add_executable(${fname} ${ldsrc})
  else ()
    cuda_add_executable(${fname} ${ldsrc})
  endif ()

  if(VTK_LIBRARIES)
    target_link_libraries(${fname} edda ${LINK_LIBS} ${VTK_LIBRARIES} netcdf_c++)
  else()
    target_link_libraries(${fname} edda ${LINK_LIBS} vtkHybrid vtkWidgets netcdf_c++)
  endif ()

  
endforeach (ldsrc)

