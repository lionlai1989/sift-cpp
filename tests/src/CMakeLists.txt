add_executable(test_image test_image.cpp)
target_link_libraries(
  test_image
  GTest::gtest_main
  img
)
target_include_directories(test_image PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_compile_features(test_image PRIVATE cxx_std_17)

add_executable(test_sift test_sift.cpp)
target_link_libraries(
  test_sift
  GTest::gtest_main
  img
  sift
)
target_include_directories(test_sift PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_compile_features(test_sift PRIVATE cxx_std_17)

include(GoogleTest)
gtest_discover_tests(test_image)
gtest_discover_tests(test_sift)
