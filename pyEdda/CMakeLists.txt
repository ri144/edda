cmake_minimum_required(VERSION 2.8.12)
project(pyEdda)

include_directories(${CMAKE_SOURCE_DIR}/src)



add_subdirectory(pybind11)

pybind11_add_module(pyedda pyEdda.cpp)

target_link_libraries(pyedda PRIVATE edda) #needed when using Visual Studio
