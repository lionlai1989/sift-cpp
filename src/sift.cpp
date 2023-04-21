#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

#include "image.hpp"
#include "sift.hpp"

namespace sift {

ScaleSpacePyramid generate_gaussian_pyramid(const Image &img, double sigma_min,
                                            int num_octaves,
                                            int scales_per_octave) {
  // assume initial sigma is 1.0 (after resizing) and smooth
  // the image with sigma_diff to reach requried base_sigma
  double base_sigma = sigma_min / MIN_PIX_DIST;
  Image base_img =
      img.resize(img.width * 2, img.height * 2, Interpolation::BILINEAR);
  double sigma_diff = std::sqrt(base_sigma * base_sigma - 1.0f);
  base_img = gaussian_blur(base_img, sigma_diff);

  int imgs_per_octave = scales_per_octave + 3;

  // determine sigma values for bluring
  double k = std::pow(2, 1.0 / scales_per_octave);
  std::vector<double> sigma_vals{base_sigma};
  for (int i = 1; i < imgs_per_octave; i++) {
    double sigma_prev = base_sigma * std::pow(k, i - 1);
    double sigma_total = k * sigma_prev;
    sigma_vals.push_back(
        std::sqrt(sigma_total * sigma_total - sigma_prev * sigma_prev));
  }

  // create a scale space pyramid of gaussian images
  // images in each octave are half the size of images in the previous one
  ScaleSpacePyramid pyramid = {num_octaves, imgs_per_octave,
                               std::vector<std::vector<Image>>(num_octaves)};
  for (int i = 0; i < num_octaves; i++) {
    pyramid.octaves[i].reserve(imgs_per_octave);
    pyramid.octaves[i].push_back(std::move(base_img));
    for (int j = 1; j < sigma_vals.size(); j++) {
      const Image &prev_img = pyramid.octaves[i].back();
      pyramid.octaves[i].push_back(gaussian_blur(prev_img, sigma_vals[j]));
    }
    // prepare base image for next octave
    const Image &next_base_img = pyramid.octaves[i][imgs_per_octave - 3];
    base_img =
        next_base_img.resize(next_base_img.width / 2, next_base_img.height / 2,
                             Interpolation::NEAREST);
  }
  return pyramid;
}

// generate pyramid of difference of gaussians (DoG) images
ScaleSpacePyramid generate_dog_pyramid(const ScaleSpacePyramid &img_pyramid) {
  ScaleSpacePyramid dog_pyramid = {
      img_pyramid.num_octaves, img_pyramid.imgs_per_octave - 1,
      std::vector<std::vector<Image>>(img_pyramid.num_octaves)};
  for (int i = 0; i < dog_pyramid.num_octaves; i++) {
    dog_pyramid.octaves[i].reserve(dog_pyramid.imgs_per_octave);
    for (int j = 1; j < img_pyramid.imgs_per_octave; j++) {
      Image diff = img_pyramid.octaves[i][j];
      diff.pixels -= img_pyramid.octaves[i][j - 1].pixels;
      dog_pyramid.octaves[i].push_back(diff);
    }
  }
  return dog_pyramid;
}

bool point_is_extremum(const std::vector<Image> &octave, int scale, int x,
                       int y) {
  const Image &img = octave[scale];
  const Image &prev = octave[scale - 1];
  const Image &next = octave[scale + 1];

  bool is_min = true, is_max = true;
  double val = img.get_pixel(x, y, 0), neighbor;

  for (int dx : {-1, 0, 1}) {
    for (int dy : {-1, 0, 1}) {
      neighbor = prev.get_pixel(x + dx, y + dy, 0);
      if (neighbor > val)
        is_max = false;
      if (neighbor < val)
        is_min = false;

      neighbor = next.get_pixel(x + dx, y + dy, 0);
      if (neighbor > val)
        is_max = false;
      if (neighbor < val)
        is_min = false;

      neighbor = img.get_pixel(x + dx, y + dy, 0);
      if (neighbor > val)
        is_max = false;
      if (neighbor < val)
        is_min = false;

      if (!is_min && !is_max)
        return false;
    }
  }
  return true;
}

// fit a quadratic near the discrete extremum,
// update the keypoint (interpolated) extremum value
// and return offsets of the interpolated extremum from the discrete extremum
std::tuple<double, double, double>
fit_quadratic(Keypoint &kp, const std::vector<Image> &octave, int scale) {
  const Image &img = octave[scale];
  const Image &prev = octave[scale - 1];
  const Image &next = octave[scale + 1];

  double g1, g2, g3;
  double h11, h12, h13, h22, h23, h33;
  int x = kp.i, y = kp.j;

  // gradient
  g1 = (next.get_pixel(x, y, 0) - prev.get_pixel(x, y, 0)) * 0.5;
  g2 = (img.get_pixel(x + 1, y, 0) - img.get_pixel(x - 1, y, 0)) * 0.5;
  g3 = (img.get_pixel(x, y + 1, 0) - img.get_pixel(x, y - 1, 0)) * 0.5;

  // hessian
  h11 = next.get_pixel(x, y, 0) + prev.get_pixel(x, y, 0) -
        2 * img.get_pixel(x, y, 0);
  h22 = img.get_pixel(x + 1, y, 0) + img.get_pixel(x - 1, y, 0) -
        2 * img.get_pixel(x, y, 0);
  h33 = img.get_pixel(x, y + 1, 0) + img.get_pixel(x, y - 1, 0) -
        2 * img.get_pixel(x, y, 0);
  h12 = (next.get_pixel(x + 1, y, 0) - next.get_pixel(x - 1, y, 0) -
         prev.get_pixel(x + 1, y, 0) + prev.get_pixel(x - 1, y, 0)) *
        0.25;
  h13 = (next.get_pixel(x, y + 1, 0) - next.get_pixel(x, y - 1, 0) -
         prev.get_pixel(x, y + 1, 0) + prev.get_pixel(x, y - 1, 0)) *
        0.25;
  h23 = (img.get_pixel(x + 1, y + 1, 0) - img.get_pixel(x + 1, y - 1, 0) -
         img.get_pixel(x - 1, y + 1, 0) + img.get_pixel(x - 1, y - 1, 0)) *
        0.25;

  // invert hessian
  double hinv11, hinv12, hinv13, hinv22, hinv23, hinv33;
  double det = h11 * h22 * h33 - h11 * h23 * h23 - h12 * h12 * h33 +
               2 * h12 * h13 * h23 - h13 * h13 * h22;
  hinv11 = (h22 * h33 - h23 * h23) / det;
  hinv12 = (h13 * h23 - h12 * h33) / det;
  hinv13 = (h12 * h23 - h13 * h22) / det;
  hinv22 = (h11 * h33 - h13 * h13) / det;
  hinv23 = (h12 * h13 - h11 * h23) / det;
  hinv33 = (h11 * h22 - h12 * h12) / det;

  // find offsets of the interpolated extremum from the discrete extremum
  double offset_s = -hinv11 * g1 - hinv12 * g2 - hinv13 * g3;
  double offset_x = -hinv12 * g1 - hinv22 * g2 - hinv23 * g3;
  double offset_y = -hinv13 * g1 - hinv23 * g3 - hinv33 * g3;

  double interpolated_extrema_val =
      img.get_pixel(x, y, 0) +
      0.5 * (g1 * offset_s + g2 * offset_x + g3 * offset_y);
  kp.extremum_val = interpolated_extrema_val;
  return {offset_s, offset_x, offset_y};
}

bool point_is_on_edge(const Keypoint &kp, const std::vector<Image> &octave,
                      double edge_thresh = C_EDGE) {
  const Image &img = octave[kp.scale];
  double h11, h12, h22;
  int x = kp.i, y = kp.j;
  h11 = img.get_pixel(x + 1, y, 0) + img.get_pixel(x - 1, y, 0) -
        2 * img.get_pixel(x, y, 0);
  h22 = img.get_pixel(x, y + 1, 0) + img.get_pixel(x, y - 1, 0) -
        2 * img.get_pixel(x, y, 0);
  h12 = (img.get_pixel(x + 1, y + 1, 0) - img.get_pixel(x + 1, y - 1, 0) -
         img.get_pixel(x - 1, y + 1, 0) + img.get_pixel(x - 1, y - 1, 0)) *
        0.25;

  double det_hessian = h11 * h22 - h12 * h12;
  double tr_hessian = h11 + h22;
  double edgeness = tr_hessian * tr_hessian / det_hessian;

  if (edgeness > std::pow(edge_thresh + 1, 2) / edge_thresh)
    return true;
  else
    return false;
}

void find_input_img_coords(Keypoint &kp, double offset_s, double offset_x,
                           double offset_y, double sigma_min = SIGMA_MIN,
                           double min_pix_dist = MIN_PIX_DIST,
                           int n_spo = N_SPO) {
  kp.sigma = std::pow(2, kp.octave) * sigma_min *
             std::pow(2, (offset_s + kp.scale) / n_spo);
  kp.x = min_pix_dist * std::pow(2, kp.octave) * (offset_x + kp.i);
  kp.y = min_pix_dist * std::pow(2, kp.octave) * (offset_y + kp.j);
}

bool refine_or_discard_keypoint(Keypoint &kp, const std::vector<Image> &octave,
                                double contrast_thresh, double edge_thresh) {
  int k = 0;
  bool kp_is_valid = false;
  while (k++ < MAX_REFINEMENT_ITERS) {
    auto [offset_s, offset_x, offset_y] = fit_quadratic(kp, octave, kp.scale);

    double max_offset =
        std::max({std::abs(offset_s), std::abs(offset_x), std::abs(offset_y)});
    // find nearest discrete coordinates
    kp.scale += std::round(offset_s);
    kp.i += std::round(offset_x);
    kp.j += std::round(offset_y);
    if (kp.scale >= octave.size() - 1 || kp.scale < 1)
      break;

    bool valid_contrast = std::abs(kp.extremum_val) > contrast_thresh;
    if (max_offset < 0.6 && valid_contrast &&
        !point_is_on_edge(kp, octave, edge_thresh)) {
      find_input_img_coords(kp, offset_s, offset_x, offset_y);
      kp_is_valid = true;
      break;
    }
  }
  return kp_is_valid;
}

std::vector<Keypoint> find_keypoints(const ScaleSpacePyramid &dog_pyramid,
                                     double contrast_thresh,
                                     double edge_thresh) {
  std::vector<Keypoint> keypoints;
  for (int i = 0; i < dog_pyramid.num_octaves; i++) {
    const std::vector<Image> &octave = dog_pyramid.octaves[i];
    for (int j = 1; j < dog_pyramid.imgs_per_octave - 1; j++) {
      const Image &img = octave[j];
      for (int x = 1; x < img.width - 1; x++) {
        for (int y = 1; y < img.height - 1; y++) {
          if (std::abs(img.get_pixel(x, y, 0)) < 0.8 * contrast_thresh) {
            continue;
          }
          if (point_is_extremum(octave, j, x, y)) {
            Keypoint kp = {x, y, i, j, -1, -1, -1, -1};
            bool kp_is_valid = refine_or_discard_keypoint(
                kp, octave, contrast_thresh, edge_thresh);
            if (kp_is_valid) {
              keypoints.push_back(kp);
            }
          }
        }
      }
    }
  }
  return keypoints;
}

// calculate x and y derivatives for all images in the input pyramid
ScaleSpacePyramid generate_gradient_pyramid(const ScaleSpacePyramid &pyramid) {
  ScaleSpacePyramid grad_pyramid = {
      pyramid.num_octaves, pyramid.imgs_per_octave,
      std::vector<std::vector<Image>>(pyramid.num_octaves)};
  for (int i = 0; i < pyramid.num_octaves; i++) {
    grad_pyramid.octaves[i].reserve(grad_pyramid.imgs_per_octave);
    int width = pyramid.octaves[i][0].width;
    int height = pyramid.octaves[i][0].height;
    for (int j = 0; j < pyramid.imgs_per_octave; j++) {
      Image grad(width, height, 2);
      double gx, gy;
      for (int x = 1; x < grad.width - 1; x++) {
        for (int y = 1; y < grad.height - 1; y++) {
          gx = (pyramid.octaves[i][j].get_pixel(x + 1, y, 0) -
                pyramid.octaves[i][j].get_pixel(x - 1, y, 0)) *
               0.5;
          grad.set_pixel(x, y, 0, gx);
          grad.pixels(0, y, x) = gx;
          gy = (pyramid.octaves[i][j].get_pixel(x, y + 1, 0) -
                pyramid.octaves[i][j].get_pixel(x, y - 1, 0)) *
               0.5;
          grad.set_pixel(x, y, 1, gy);
          grad.pixels(1, y, x) = gy;
        }
      }
      grad_pyramid.octaves[i].push_back(grad);
    }
  }
  return grad_pyramid;
}

// convolve 6x with box filter
void smooth_histogram(double hist[N_BINS]) {
  double tmp_hist[N_BINS];
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < N_BINS; j++) {
      int prev_idx = (j - 1 + N_BINS) % N_BINS;
      int next_idx = (j + 1) % N_BINS;
      tmp_hist[j] = (hist[prev_idx] + hist[j] + hist[next_idx]) / 3;
    }
    for (int j = 0; j < N_BINS; j++) {
      hist[j] = tmp_hist[j];
    }
  }
}

std::vector<double>
find_keypoint_orientations(Keypoint &kp, const ScaleSpacePyramid &grad_pyramid,
                           double lambda_ori, double lambda_desc) {
  double pix_dist = MIN_PIX_DIST * std::pow(2, kp.octave);
  const Image &img_grad = grad_pyramid.octaves[kp.octave][kp.scale];

  // discard kp if too close to image borders
  double min_dist_from_border =
      std::min({kp.x, kp.y, pix_dist * img_grad.width - kp.x,
                pix_dist * img_grad.height - kp.y});
  if (min_dist_from_border <= std::sqrt(2) * lambda_desc * kp.sigma) {
    return {};
  }

  double hist[N_BINS] = {};
  int bin;
  double gx, gy, grad_norm, weight, theta;
  double patch_sigma = lambda_ori * kp.sigma;
  double patch_radius = 3 * patch_sigma;
  int x_start = std::round((kp.x - patch_radius) / pix_dist);
  int x_end = std::round((kp.x + patch_radius) / pix_dist);
  int y_start = std::round((kp.y - patch_radius) / pix_dist);
  int y_end = std::round((kp.y + patch_radius) / pix_dist);

  // accumulate gradients in orientation histogram
  for (int x = x_start; x <= x_end; x++) {
    for (int y = y_start; y <= y_end; y++) {
      gx = img_grad.get_pixel(x, y, 0);
      gy = img_grad.get_pixel(x, y, 1);
      grad_norm = std::sqrt(gx * gx + gy * gy);
      weight = std::exp(-(std::pow(x * pix_dist - kp.x, 2) +
                          std::pow(y * pix_dist - kp.y, 2)) /
                        (2 * patch_sigma * patch_sigma));
      theta = std::fmod(std::atan2(gy, gx) + 2 * M_PI, 2 * M_PI);
      bin = (int)std::round(N_BINS / (2 * M_PI) * theta) % N_BINS;
      hist[bin] += weight * grad_norm;
    }
  }

  smooth_histogram(hist);

  // extract reference orientations
  double ori_thresh = 0.8, ori_max = 0;
  std::vector<double> orientations;
  for (int j = 0; j < N_BINS; j++) {
    if (hist[j] > ori_max) {
      ori_max = hist[j];
    }
  }
  for (int j = 0; j < N_BINS; j++) {
    if (hist[j] >= ori_thresh * ori_max) {
      double prev = hist[(j - 1 + N_BINS) % N_BINS],
             next = hist[(j + 1) % N_BINS];
      if (prev > hist[j] || next > hist[j])
        continue;
      double theta =
          2 * M_PI * (j + 1) / N_BINS +
          M_PI / N_BINS * (prev - next) / (prev - 2 * hist[j] + next);
      orientations.push_back(theta);
    }
  }
  return orientations;
}

void update_histograms(double hist[N_HIST][N_HIST][N_ORI], double x, double y,
                       double contrib, double theta_mn, double lambda_desc) {
  double x_i, y_j;
  for (int i = 1; i <= N_HIST; i++) {
    x_i = (i - (1 + (double)N_HIST) / 2) * 2 * lambda_desc / N_HIST;
    if (std::abs(x_i - x) > 2 * lambda_desc / N_HIST)
      continue;
    for (int j = 1; j <= N_HIST; j++) {
      y_j = (j - (1 + (double)N_HIST) / 2) * 2 * lambda_desc / N_HIST;
      if (std::abs(y_j - y) > 2 * lambda_desc / N_HIST)
        continue;

      double hist_weight =
          (1 - N_HIST * 0.5 / lambda_desc * std::abs(x_i - x)) *
          (1 - N_HIST * 0.5 / lambda_desc * std::abs(y_j - y));

      for (int k = 1; k <= N_ORI; k++) {
        double theta_k = 2 * M_PI * (k - 1) / N_ORI;
        double theta_diff = std::fmod(theta_k - theta_mn + 2 * M_PI, 2 * M_PI);
        if (std::abs(theta_diff) >= 2 * M_PI / N_ORI)
          continue;
        double bin_weight = 1 - N_ORI * 0.5 / M_PI * std::abs(theta_diff);
        hist[i - 1][j - 1][k - 1] += hist_weight * bin_weight * contrib;
      }
    }
  }
}

void hists_to_vec(double histograms[N_HIST][N_HIST][N_ORI],
                  std::array<uint8_t, 128> &feature_vec) {
  int size = N_HIST * N_HIST * N_ORI;
  double *hist = reinterpret_cast<double *>(histograms);

  double norm = 0;
  for (int i = 0; i < size; i++) {
    norm += hist[i] * hist[i];
  }
  norm = std::sqrt(norm);
  double norm2 = 0;
  for (int i = 0; i < size; i++) {
    hist[i] = std::min(hist[i], 0.2f * norm);
    norm2 += hist[i] * hist[i];
  }
  norm2 = std::sqrt(norm2);
  for (int i = 0; i < size; i++) {
    double val = std::floor(512 * hist[i] / norm2);
    feature_vec[i] = std::min((int)val, 255);
  }
}

void compute_keypoint_descriptor(Keypoint &kp, double theta,
                                 const ScaleSpacePyramid &grad_pyramid,
                                 double lambda_desc) {
  double pix_dist = MIN_PIX_DIST * std::pow(2, kp.octave);
  const Image &img_grad = grad_pyramid.octaves[kp.octave][kp.scale];
  double histograms[N_HIST][N_HIST][N_ORI] = {0};

  // find start and end coords for loops over image patch
  double half_size =
      std::sqrt(2) * lambda_desc * kp.sigma * (N_HIST + 1.) / N_HIST;
  int x_start = std::round((kp.x - half_size) / pix_dist);
  int x_end = std::round((kp.x + half_size) / pix_dist);
  int y_start = std::round((kp.y - half_size) / pix_dist);
  int y_end = std::round((kp.y + half_size) / pix_dist);

  double cos_t = std::cos(theta), sin_t = std::sin(theta);
  double patch_sigma = lambda_desc * kp.sigma;
  // accumulate samples into histograms
  for (int m = x_start; m <= x_end; m++) {
    for (int n = y_start; n <= y_end; n++) {
      // find normalized coords w.r.t. kp position and reference orientation
      double x =
          ((m * pix_dist - kp.x) * cos_t + (n * pix_dist - kp.y) * sin_t) /
          kp.sigma;
      double y =
          (-(m * pix_dist - kp.x) * sin_t + (n * pix_dist - kp.y) * cos_t) /
          kp.sigma;

      // verify (x, y) is inside the description patch
      if (std::max(std::abs(x), std::abs(y)) >
          lambda_desc * (N_HIST + 1.) / N_HIST)
        continue;

      double gx = img_grad.get_pixel(m, n, 0), gy = img_grad.get_pixel(m, n, 1);
      double theta_mn =
          std::fmod(std::atan2(gy, gx) - theta + 4 * M_PI, 2 * M_PI);
      double grad_norm = std::sqrt(gx * gx + gy * gy);
      double weight = std::exp(-(std::pow(m * pix_dist - kp.x, 2) +
                                 std::pow(n * pix_dist - kp.y, 2)) /
                               (2 * patch_sigma * patch_sigma));
      double contribution = weight * grad_norm;

      update_histograms(histograms, x, y, contribution, theta_mn, lambda_desc);
    }
  }

  // build feature vector (descriptor) from histograms
  hists_to_vec(histograms, kp.descriptor);
}

std::vector<Keypoint>
find_keypoints_and_descriptors(const Image &img, double sigma_min,
                               int num_octaves, int scales_per_octave,
                               double contrast_thresh, double edge_thresh,
                               double lambda_ori, double lambda_desc) {
  assert(img.channels == 1 || img.channels == 3);

  const Image &input = img.channels == 1 ? img : rgb_to_grayscale(img);
  ScaleSpacePyramid gaussian_pyramid = generate_gaussian_pyramid(
      input, sigma_min, num_octaves, scales_per_octave);
  ScaleSpacePyramid dog_pyramid = generate_dog_pyramid(gaussian_pyramid);
  std::vector<Keypoint> tmp_kps =
      find_keypoints(dog_pyramid, contrast_thresh, edge_thresh);
  ScaleSpacePyramid grad_pyramid = generate_gradient_pyramid(gaussian_pyramid);

  std::vector<Keypoint> kps;

  for (Keypoint &kp_tmp : tmp_kps) {
    std::vector<double> orientations = find_keypoint_orientations(
        kp_tmp, grad_pyramid, lambda_ori, lambda_desc);
    for (double theta : orientations) {
      Keypoint kp = kp_tmp;
      compute_keypoint_descriptor(kp, theta, grad_pyramid, lambda_desc);
      kps.push_back(kp);
    }
  }
  return kps;
}

double euclidean_dist(std::array<uint8_t, 128> &a,
                      std::array<uint8_t, 128> &b) {
  double dist = 0;
  for (int i = 0; i < 128; i++) {
    int di = (int)a[i] - b[i];
    dist += di * di;
  }
  return std::sqrt(dist);
}

std::vector<std::pair<int, int>> find_keypoint_matches(std::vector<Keypoint> &a,
                                                       std::vector<Keypoint> &b,
                                                       double thresh_relative,
                                                       double thresh_absolute) {
  assert(a.size() >= 2 && b.size() >= 2);

  std::vector<std::pair<int, int>> matches;

  for (int i = 0; i < a.size(); i++) {
    // find two nearest neighbours in b for current keypoint from a
    int nn1_idx = -1;
    double nn1_dist = 100000000, nn2_dist = 100000000;
    for (int j = 0; j < b.size(); j++) {
      double dist = euclidean_dist(a[i].descriptor, b[j].descriptor);
      if (dist < nn1_dist) {
        nn2_dist = nn1_dist;
        nn1_dist = dist;
        nn1_idx = j;
      } else if (nn1_dist <= dist && dist < nn2_dist) {
        nn2_dist = dist;
      }
    }
    if (nn1_dist < thresh_relative * nn2_dist && nn1_dist < thresh_absolute) {
      matches.push_back({i, nn1_idx});
    }
  }
  return matches;
}

Image draw_keypoints(const Image &img, const std::vector<Keypoint> &kps) {
  Image res(img);
  if (img.channels == 1) {
    res = grayscale_to_rgb(res);
  }
  for (auto &kp : kps) {
    draw_point(res, kp.x, kp.y, 5);
  }
  return res;
}

Image draw_matches(const Image &a, const Image &b, std::vector<Keypoint> &kps_a,
                   std::vector<Keypoint> &kps_b,
                   std::vector<std::pair<int, int>> matches) {
  Image res(a.width + b.width, std::max(a.height, b.height), 3);

  for (int i = 0; i < a.width; i++) {
    for (int j = 0; j < a.height; j++) {
      res.set_pixel(i, j, 0, a.get_pixel(i, j, 0));
      res.set_pixel(i, j, 1, a.get_pixel(i, j, a.channels == 3 ? 1 : 0));
      res.set_pixel(i, j, 2, a.get_pixel(i, j, a.channels == 3 ? 2 : 0));
      res.pixels(0, j, i) = a.get_pixel(i, j, 0);
      res.pixels(1, j, i) = a.get_pixel(i, j, a.channels == 3 ? 1 : 0);
      res.pixels(2, j, i) = a.get_pixel(i, j, a.channels == 3 ? 2 : 0);
    }
  }
  for (int i = 0; i < b.width; i++) {
    for (int j = 0; j < b.height; j++) {
      res.set_pixel(a.width + i, j, 0, b.get_pixel(i, j, 0));
      res.set_pixel(a.width + i, j, 1,
                    b.get_pixel(i, j, b.channels == 3 ? 1 : 0));
      res.set_pixel(a.width + i, j, 2,
                    b.get_pixel(i, j, b.channels == 3 ? 2 : 0));
      res.pixels(0, j, a.width + i) = b.get_pixel(i, j, 0);
      res.pixels(1, j, a.width + i) =
          b.get_pixel(i, j, b.channels == 3 ? 1 : 0);
      res.pixels(2, j, a.width + i) =
          b.get_pixel(i, j, b.channels == 3 ? 2 : 0);
    }
  }

  for (auto &m : matches) {
    Keypoint &kp_a = kps_a[m.first];
    Keypoint &kp_b = kps_b[m.second];
    draw_line(res, kp_a.x, kp_a.y, a.width + kp_b.x, kp_b.y);
  }
  return res;
}

} // namespace sift
