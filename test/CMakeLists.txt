project (test)


add_executable (testGaussian testGaussian.C )
target_link_libraries(testGaussian edda)

add_executable (testVectorField testVectorField.C )
target_link_libraries(testVectorField edda)

