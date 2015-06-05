# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

project (edda_test)


set(TEST_SOURCES
    test_gaussian.cpp
    )

foreach (test ${TEST_SOURCES})
  get_filename_component(tname ${test} NAME_WE)
  add_executable(${tname} ${test})
  target_link_libraries(${tname} edda)

  add_test(NAME ${tname}
    COMMAND ${tname}
    )
endforeach (test)
