#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>

#include <utility>
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include "image.hpp"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

Image::Image(const xt::xtensor<double, 3> &input_matrix) {
  // input_matrix: channel, height, width
  this->width = input_matrix.shape(2);
  this->height = input_matrix.shape(1);
  this->channels = input_matrix.shape(0);
  pixels = xt::xtensor<double, 3>(input_matrix);

  // auto max_value = input_matrix.maxCoeff();

  this->size = this->width * this->height * this->channels;
}

Image::Image(std::string file_path) {
  using ImageDate =
      std::unique_ptr<unsigned char[], decltype(&stbi_image_free)>;
  ImageDate img_data{
      stbi_load(file_path.c_str(), &width, &height, &channels, 0),
      stbi_image_free};

  if (img_data == nullptr) {
    const char *error_msg = stbi_failure_reason();
    std::string err_msg = "Failed to load image: " + file_path + '\n' +
                          "Error msg (stb_image): " + error_msg + '\n';
    throw std::runtime_error(err_msg.c_str());
  }

  pixels = xt::zeros<double>({channels, height, width});

  this->size = this->width * this->height * this->channels;

  for (int x = 0; x < this->width; x++) {
    for (int y = 0; y < this->height; y++) {
      for (int c = 0; c < this->channels; c++) {
        int src_idx = y * this->width * this->channels + x * this->channels + c;
        pixels(c, y, x) = img_data[src_idx] / 255.;
      }
    }
  }
  if (this->channels == 4) {
    this->channels = 3; // ignore alpha channel
  }
}

Image::Image(int w, int h, int c)
    : width{w}, height{h}, channels{c}, size{w * h * c}, pixels{
                                                             xt::zeros<double>(
                                                                 {c, h, w})} {}

Image::Image() : width{0}, height{0}, channels{0}, size{0} {}

Image::~Image() {}

Image::Image(const Image &other)
    : width{other.width}, height{other.height}, channels{other.channels},
      size{other.size}, pixels{xt::zeros<double>(
                            {other.channels, other.height, other.width})} {
  // std::cout << "copy constructor\n";
  pixels = other.pixels;
  //   for (int i = 0; i < size; i++)
  //     data[i] = other.data[i];
}

Image::Image(Image &&other)
    : width{other.width}, height{other.height}, channels{other.channels},
      size{other.size}, pixels{std::move(other.pixels)} {
  // std::cout << "move constructor\n";
}

Image &Image::operator=(Image rhs) {
  _swap(rhs);
  return *this;
}

void Image::_swap(Image &other) {
  std::swap(channels, other.channels);
  std::swap(height, other.height);
  std::swap(width, other.width);
  std::swap(size, other.size);
  std::swap(pixels, other.pixels);
}

// save image as jpg file
bool Image::save(std::string file_path) {
  auto file_extension = std::filesystem::path(file_path).extension();

  std::unique_ptr<unsigned char[]> out_data =
      std::make_unique<unsigned char[]>(width * height * channels);

  for (auto x = 0; x < width; ++x) {
    for (auto y = 0; y < height; ++y) {
      for (auto c = 0; c < channels; ++c) {
        int dst_idx = y * width * channels + x * channels + c;
        // Fill out_data with uint8 values range 0-255
        out_data[dst_idx] = std::lroundf(pixels(c, y, x) * 255.0);
      }
    }
  }

  bool success{false};
  if (file_extension == std::string(".jpg") ||
      file_extension == std::string(".JPG")) {
    auto quality = 100;
    success = stbi_write_jpg(file_path.c_str(), width, height, channels,
                             out_data.get(), quality);
  } else if (file_extension == std::string(".png") ||
             file_extension == std::string(".png")) {
    auto stride_in_bytes = width * channels;
    success = stbi_write_png(file_path.c_str(), width, height, channels,
                             out_data.get(), stride_in_bytes);
  } else {
    std::cerr << "Unsupported file format: " << file_extension << "\n";
  }
  if (!success)
    std::cerr << "Failed to save image: " << file_path << "\n";

  return true;
}

void Image::set_pixel(int x, int y, int c, double val) {
  if (x >= this->width || x < 0 || y >= this->height || y < 0 ||
      c >= this->channels || c < 0) {
    std::cerr << "set_pixel() error: Index out of bounds.\n";
    std::exit(1);
  }
}

double Image::get_pixel(int x, int y, int c) const {
  if (x < 0)
    x = 0;
  if (x >= this->width)
    x = this->width - 1;
  if (y < 0)
    y = 0;
  if (y >= this->height)
    y = this->height - 1;
  return pixels(c, y, x);
}

void Image::clamp() {
  int size = this->width * this->height * this->channels;
  for (int i = 0; i < size; i++) {
    double val = 0;
    // double val = this->data[i];
    val = (val > 1.0) ? 1.0 : val;
    val = (val < 0.0) ? 0.0 : val;
    // this->data[i] = val;
  }
}

bool Image::operator==(const Image &other) const {
  return (width == other.width) && (height == other.height) &&
         (channels == other.channels) && (size == other.size) &&
         (pixels == other.pixels);
}

// map coordinate from 0-current_max range to 0-new_max range
double myround(double d) { return std::floor(d + 0.5); }
double map_coordinate(double new_max, double current_max, double coord) {
  //   double a = new_max / current_max;
  //   double b = -0.5 + a * 0.5;
  //   return a * coord + b;
  double slope = 1.0 * (new_max) / (current_max);
  return myround(slope * (coord));
}

Image Image::resize(int new_w, int new_h, Interpolation method) const {
  Image resized(new_w, new_h, this->channels);
  double value = 0;
  for (int x = 0; x < new_w; x++) {
    for (int y = 0; y < new_h; y++) {
      for (int c = 0; c < resized.channels; c++) {
        assert(x >= 0 && x < new_w && y >= 0 && y < new_h);
        // NOTE: Array with 100 elements, its max index is 100-1.
        double old_x = map_coordinate(this->width - 1, new_w - 1, x);
        double old_y = map_coordinate(this->height - 1, new_h - 1, y);
        if (method == Interpolation::BILINEAR) {
          value = bilinear_interpolate(*this, old_x, old_y, c);
        } else if (method == Interpolation::NEAREST) {
          value = nn_interpolate(*this, old_x, old_y, c);
        }
        resized.set_pixel(x, y, c, value);
        resized.pixels(c, y, x) = value;
      }
    }
  }
  return resized;
}

double bilinear_interpolate(const Image &img, double x, double y, int c) {
  double p1, p2, p3, p4, q1, q2;
  double x_floor = std::floor(x);
  double y_floor = std::floor(y);

  double x_ceil = (x_floor + 1) > (img.width - 1) ? x_floor : x_floor + 1;
  double y_ceil = (y_floor + 1) > (img.height - 1) ? y_floor : y_floor + 1;
  assert(x_floor >= 0 && x_floor < img.width && x_ceil >= 0 &&
         x_ceil < img.width && y_floor >= 0 && y_floor < img.height &&
         y_ceil >= 0 && y_ceil < img.height);
  p1 = img.get_pixel(x_floor, y_floor, c);
  p2 = img.get_pixel(x_ceil, y_floor, c);
  p3 = img.get_pixel(x_floor, y_ceil, c);
  p4 = img.get_pixel(x_ceil, y_ceil, c);
  q1 = (y_ceil - y) * p1 + (y - y_floor) * p3;
  q2 = (y_ceil - y) * p2 + (y - y_floor) * p4;
  return (x_ceil - x) * q1 + (x - x_floor) * q2;
}

double nn_interpolate(const Image &img, double x, double y, int c) {

  int round_x =
      std::lround(x) > (img.width - 1) ? (img.width - 1) : std::lround(x);
  int round_y =
      std::lround(y) > (img.height - 1) ? (img.height - 1) : std::lround(y);
  if (round_x < 0) {
    round_x = 0;
  }
  if (round_y < 0) {
    round_y = 0;
  }
  assert(round_x >= 0);
  assert(round_x < img.width);
  assert(round_y >= 0);
  assert(round_y < img.height);
  return img.get_pixel(round_x, round_y, c);
}

Image rgb_to_grayscale(const Image &img) {
  assert(img.channels == 3);
  Image gray(img.width, img.height, 1);
  for (int x = 0; x < img.width; x++) {
    for (int y = 0; y < img.height; y++) {
      double red, green, blue;
      red = img.get_pixel(x, y, 0);
      green = img.get_pixel(x, y, 1);
      blue = img.get_pixel(x, y, 2);
      gray.set_pixel(x, y, 0, 0.299 * red + 0.587 * green + 0.114 * blue);
    }
  }

  xt::xtensor<double, 2> red = xt::view(img.pixels, 0, xt::all(), xt::all());
  xt::xtensor<double, 2> green = xt::view(img.pixels, 1, xt::all(), xt::all());
  xt::xtensor<double, 2> blue = xt::view(img.pixels, 2, xt::all(), xt::all());
  xt::view(gray.pixels, 0, xt::all(), xt::all()) =
      0.299 * red + 0.587 * green + 0.114 * blue;
  return gray;
}

Image grayscale_to_rgb(const Image &img) {
  assert(img.channels == 1);
  Image rgb(img.width, img.height, 3);
  for (int x = 0; x < img.width; x++) {
    for (int y = 0; y < img.height; y++) {
      double gray_val = img.get_pixel(x, y, 0);
      rgb.set_pixel(x, y, 0, gray_val);
      rgb.set_pixel(x, y, 1, gray_val);
      rgb.set_pixel(x, y, 2, gray_val);
      rgb.pixels(0, y, x) = gray_val;
      rgb.pixels(1, y, x) = gray_val;
      rgb.pixels(2, y, x) = gray_val;
    }
  }
  return rgb;
}

// separable 2D gaussian blur for 1 channel image
Image gaussian_blur(const Image &img, double sigma) {
  assert(img.channels == 1);

  int size = std::ceil(6 * sigma);
  if (size % 2 == 0)
    size++;
  int center = size / 2;
  Image kernel(size, 1, 1);
  double sum = 0;
  for (int k = -size / 2; k <= size / 2; k++) {
    double val = std::exp(-(k * k) / (2 * sigma * sigma));
    assert(center + k >= 0 && center + k < size);
    kernel.set_pixel(center + k, 0, 0, val);
    kernel.pixels(0, 0, center + k) = val;
    sum += val;
  }
  for (int k = 0; k < size; k++) {
    // kernel.data[k] /= sum;
    kernel.pixels(0, 0, k) /= sum;
  }

  Image tmp(img.width, img.height, 1);
  Image filtered(img.width, img.height, 1);

  // convolve vertical
  for (int x = 0; x < img.width; x++) {
    for (int y = center; y < img.height - center; y++) {
      double sum = 0;
      for (int k = 0; k < size; k++) {
        int dy = -center + k;
        assert((y + dy) >= 0 && (y + dy) < img.height);
        // sum += (img.get_pixel(x, y + dy, 0) * kernel.data[k]);
        sum += (img.get_pixel(x, y + dy, 0) * kernel.pixels(0, 0, k));
      }
      tmp.set_pixel(x, y, 0, sum);
      tmp.pixels(0, y, x) = sum;
    }
  }
  // convolve horizontal
  for (int x = center; x < img.width - center; x++) {
    for (int y = 0; y < img.height; y++) {
      double sum = 0;
      for (int k = 0; k < size; k++) {
        int dx = -center + k;
        assert((x + dx) >= 0 && (x + dx) < img.width);
        // sum += (tmp.get_pixel(x + dx, y, 0) * kernel.data[k]);
        sum += (tmp.get_pixel(x + dx, y, 0) * kernel.pixels(0, 0, k));
      }
      filtered.set_pixel(x, y, 0, sum);
      filtered.pixels(0, y, x) = sum;
    }
  }
  return filtered;
}

void draw_point(Image &img, int x, int y, int size) {
  for (int i = x - size / 2; i <= x + size / 2; i++) {
    for (int j = y - size / 2; j <= y + size / 2; j++) {
      if (i < 0 || i >= img.width)
        continue;
      if (j < 0 || j >= img.height)
        continue;
      if (std::abs(i - x) + std::abs(j - y) > size / 2)
        continue;
      assert(i >= 0 && j >= 0);
      assert(i < img.width && j < img.height);
      if (img.channels == 3) {
        img.set_pixel(i, j, 0, 1.f);
        img.set_pixel(i, j, 1, 0.f);
        img.set_pixel(i, j, 2, 0.f);
        img.pixels(0, j, i) = 1.f;
        img.pixels(1, j, i) = 0.f;
        img.pixels(2, j, i) = 0.f;
      } else {
        img.set_pixel(i, j, 0, 1.f);
        img.pixels(0, j, i) = 1.f;
      }
    }
  }
}

void draw_line(Image &img, int x1, int y1, int x2, int y2) {
  if (x2 < x1) {
    std::swap(x1, x2);
    std::swap(y1, y2);
  }
  int dx = x2 - x1, dy = y2 - y1;
  for (int x = x1; x < x2; x++) {
    int y = y1 + dy * (x - x1) / dx;
    assert(x >= 0 && y >= 0);
    assert(x < img.width && y < img.height);
    if (img.channels == 3) {
      img.set_pixel(x, y, 0, 0.f);
      img.set_pixel(x, y, 1, 1.f);
      img.set_pixel(x, y, 2, 0.f);
      img.pixels(0, y, x) = 0.f;
      img.pixels(1, y, x) = 1.f;
      img.pixels(2, y, x) = 0.f;
    } else {
      img.set_pixel(x, y, 0, 1.f);
      img.pixels(0, y, x) = 1.f;
    }
  }
}
