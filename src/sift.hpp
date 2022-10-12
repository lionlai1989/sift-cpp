
#ifndef SIFT_H
#define SIFT_H

#include <vector>
#include <array>
#include <cstdint>

#include "image.hpp"

namespace sift {

struct ScaleSpacePyramid {
    int num_octaves;
    int imgs_per_octave;
    std::vector<std::vector<Image>> octaves; 
};

struct Keypoint {
    // discrete coordinates
    int i;
    int j;
    int octave;
    int scale; //index of gaussian image inside the octave

    // continuous coordinates (interpolated)
    double x;
    double y;
    double sigma;
    double extremum_val; //value of interpolated DoG extremum
    
    std::array<uint8_t, 128> descriptor;
};

//*******************************************
// SIFT algorithm parameters, used by default
//*******************************************

// digital scale space configuration and keypoint detection
const int MAX_REFINEMENT_ITERS = 5;
const double SIGMA_MIN = 0.8;
const double MIN_PIX_DIST = 0.5;
const double SIGMA_IN = 0.5;
const int N_OCT = 8;
const int N_SPO = 3;
const double C_DOG = 0.015;
const double C_EDGE = 10;

// computation of the SIFT descriptor
const int N_BINS = 36;
const double LAMBDA_ORI = 1.5;
const int N_HIST = 4;
const int N_ORI = 8;
const double LAMBDA_DESC = 6;

// feature matching
const double THRESH_ABSOLUTE = 350;
const double THRESH_RELATIVE = 0.7;

ScaleSpacePyramid generate_gaussian_pyramid(const Image& img, double sigma_min=SIGMA_MIN,
                                            int num_octaves=N_OCT, int scales_per_octave=N_SPO);

ScaleSpacePyramid generate_dog_pyramid(const ScaleSpacePyramid& img_pyramid);

std::vector<Keypoint> find_keypoints(const ScaleSpacePyramid& dog_pyramid,
                                     double contrast_thresh=C_DOG, double edge_thresh=C_EDGE);

ScaleSpacePyramid generate_gradient_pyramid(const ScaleSpacePyramid& pyramid);

std::vector<double> find_keypoint_orientations(Keypoint& kp, const ScaleSpacePyramid& grad_pyramid,
                                              double lambda_ori=LAMBDA_ORI, double lambda_desc=LAMBDA_DESC);

void compute_keypoint_descriptor(Keypoint& kp, double theta, const ScaleSpacePyramid& grad_pyramid,
                                 double lambda_desc=LAMBDA_DESC);

std::vector<Keypoint> find_keypoints_and_descriptors(const Image& img, double sigma_min=SIGMA_MIN,
                                                     int num_octaves=N_OCT, 
                                                     int scales_per_octave=N_SPO, 
                                                     double contrast_thresh=C_DOG,
                                                     double edge_thresh=C_EDGE,
                                                     double lambda_ori=LAMBDA_ORI,
                                                     double lambda_desc=LAMBDA_DESC);

std::vector<std::pair<int, int>> find_keypoint_matches(std::vector<Keypoint>& a,
                                                       std::vector<Keypoint>& b,
                                                       double thresh_relative=THRESH_RELATIVE,
                                                       double thresh_absolute=THRESH_ABSOLUTE);

Image draw_keypoints(const Image& img, const std::vector<Keypoint>& kps);

Image draw_matches(const Image& a, const Image& b, std::vector<Keypoint>& kps_a,
                   std::vector<Keypoint>& kps_b, std::vector<std::pair<int, int>> matches);

} // namespace sift
#endif
