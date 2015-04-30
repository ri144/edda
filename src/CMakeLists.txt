# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

file (GLOB SRC_FILES "*.cpp")

add_library (${PROJECT_NAME} ${LIB_TYPE} ${SRC_FILES})

set_target_properties (${PROJECT_NAME} PROPERTIES
    DEBUG_OUTPUT_NAME "${PROJECT_NAME}_d"
    RELEASE_OUTPUT_NAME "${PROJECT_NAME}")

install (
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION include/edda
    FILES_MATCHING PATTERN "*.h")

install (
    TARGETS ${PROJECT_NAME}
    DESTINATION lib/edda)
