# include_directories(../extern/eigen-3.4.0)

add_library(stb_image STATIC
            stb_image.h
            stb_image_write.h
)

set_target_properties(stb_image PROPERTIES LINKER_LANGUAGE CXX)

add_library(img STATIC
            image.hpp
            image.cpp
)
target_include_directories(img PUBLIC ${CMAKE_SOURCE_DIR}/extern)
target_link_libraries(img PUBLIC xtensor xtensor::optimize xtensor::use_xsimd stb_image)
target_compile_features(img PRIVATE cxx_std_17)

# add_library(sift STATIC
#             sift.cpp
#             sift.hpp
# )

# target_link_libraries(sift PUBLIC
#                       img
# )
