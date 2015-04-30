# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

cmake_minimum_required (VERSION 2.8)
project (edda)

file (GLOB HEADER_FILES "*.h")
file (GLOB SRC_FILES "*.C")

add_library(${PROJECT_NAME} ${LIB_TYPE} ${SRC_FILES})

set_target_properties(${PROJECT_NAME} PROPERTIES
    DEBUG_OUTPUT_NAME "${PROJECT_NAME}_d"
    RELEASE_OUTPUT_NAME "${PROJECT_NAME}"
    )

install (
    FILES ${HEADER_FILES}
    DESTINATION include/edda)

install (
    TARGETS ${PROJECT_NAME}
    DESTINATION lib/edda)


