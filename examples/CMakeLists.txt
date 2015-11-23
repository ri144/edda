# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

set(sources
	loadGaussianScalarField
	loadGaussianVectorField
    uncertainIsocontour
	streamlineGaussian
)


foreach (file ${sources})
    get_filename_component(name ${file} NAME_WE)
    add_executable ( ${name} ${file} )
    target_link_libraries ( ${name} edda )
endforeach ()

