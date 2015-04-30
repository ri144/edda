# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

add_executable(gaussian_test gaussian_test.cpp)
target_link_libraries(gaussian_test edda)

add_executable(test_vector_field test_vector_field.cpp)
target_link_libraries(test_vector_field edda)
