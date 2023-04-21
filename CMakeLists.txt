cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall --std=c++17 -O3 -fPIC" )
# set( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -O3 -fPIC" )

project(sift_cpp LANGUAGES CXX)

list(APPEND CMAKE_PREFIX_PATH "/tmp/xtl-install")
list(APPEND CMAKE_PREFIX_PATH "/tmp/xtensor-install")
list(APPEND CMAKE_PREFIX_PATH "/tmp/xsimd-install")

add_subdirectory(extern/pybind11)

# find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)

# If using `add_subdirectory` on `xtl`, `xtensor`, and `xsimd`, there will be error.
# Check https://github.com/xtensor-stack/xtensor-python/issues/293
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xsimd REQUIRED)

# add_subdirectory(extern/xtensor-python)

add_subdirectory(src)
# add_subdirectory(examples)

# include_directories(./extern/eigen-3.4.0)

# pybind11_add_module(pybind11_wrapper ./pybind11_wrapper.cpp)
# Figure out what's the difference of include_directories vs target_include_directories
# target_include_directories(pybind11_wrapper PRIVATE ./src)
# target_link_libraries(pybind11_wrapper PUBLIC img sift)

# It must be called before adding `tests` subdirectory.
enable_testing()
add_subdirectory(tests)
