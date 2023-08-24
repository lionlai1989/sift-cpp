#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "sift.hpp"
#include "image.hpp"

#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>

#define FORCE_IMPORT_ARRAY

#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pyvectorize.hpp>

namespace py = pybind11;

std::vector<sift::Keypoint> computeKeypointsAndDescriptors(const xt::xarray<double> &pixels) {
    Image input_image(pixels);
    std::vector<sift::Keypoint> kps = sift::find_keypoints_and_descriptors(input_image);
    return kps;
}

PYBIND11_MODULE(lion_sift_cpp, m) {
    m.def("computeKeypointsAndDescriptors", &computeKeypointsAndDescriptors);
    py::class_<sift::Keypoint>(m, "sift::Keypoint")
    .def(py::init<>()) // <-- bind the default constructor
    .def_readwrite("i", &sift::Keypoint::i)
    .def_readwrite("j", &sift::Keypoint::j)
    .def_readwrite("octave", &sift::Keypoint::octave)
    .def_readwrite("scale", &sift::Keypoint::scale)
    .def_readwrite("x", &sift::Keypoint::x)
    .def_readwrite("y", &sift::Keypoint::y)
    .def_readwrite("sigma", &sift::Keypoint::sigma)
    .def_readwrite("extremum_val", &sift::Keypoint::extremum_val)
    .def_readwrite("descriptor", &sift::Keypoint::descriptor);
}
