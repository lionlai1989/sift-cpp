cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall --std=c++17 -O3 -fPIC" )
# set( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -O3 -fPIC" )

project(sift_cpp LANGUAGES CXX)


add_subdirectory(src)
add_subdirectory(examples)

include_directories(./extern/eigen-3.4.0)

add_subdirectory(./extern/pybind11)
pybind11_add_module(pybind11_wrapper ./pybind11_wrapper.cpp)
target_include_directories(pybind11_wrapper PRIVATE ./src)
target_link_libraries(pybind11_wrapper PUBLIC img sift)
