# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

add_executable(test_gaussian test_gaussian.cpp)
target_link_libraries(test_gaussian edda)

add_executable(test_vector_field test_vector_field.cpp)
target_link_libraries(test_vector_field edda)

add_executable(test_gaussian_scalar_field test_gaussian_scalar_field.cpp)
target_link_libraries(test_gaussian_scalar_field edda)
