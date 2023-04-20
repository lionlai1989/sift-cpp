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

Image::Image(const xt::xtensor<double, 3> &input_matrix)
    : channels(input_matrix.shape(0)), height(input_matrix.shape(1)),
      width(input_matrix.shape(2)),
      size(input_matrix.size()), pixels{xt::xtensor<double, 3>(input_matrix)} {
  std::clog << "The Constructor takes xtensor.\n";
}

Image::Image(std::string file_path) {
  std::clog << "The constructor takes a file path.\n";

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

  size = pixels.size();
  std::clog << "The image shape: " << channels << " x " << height << " x "
            << width << '\n';
  assert(size == channels * height * width);

  for (int x = 0; x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      for (int c = 0; c < channels; ++c) {
        std::size_t src_idx = y * width * channels + x * channels + c;
        // Rescale uint8 to float 0-1.0
        pixels(c, y, x) = img_data[src_idx] / 255.0;
      }
    }
  }
  if (channels == 4) {
    channels = 3; // ignore alpha channel
  }
}

Image::Image(int c, int h, int w)
    : channels{c}, height{h}, width{w}, size{c * h * w},
      pixels{xt::xtensor<double, 3>(xt::zeros<double>({c, h, w}))} {
  std::clog << "The constructor takes (c, h, w).\n";
}

Image::Image() : channels{0}, height{0}, width{0}, size{0} {
  std::clog << "The default constructor takes no paramters.\n";
  /**
   * pixels is not initialized because its shape is unknown.
   */
}

Image::~Image() { std::clog << "Destruct Image.\n"; }

bool Image::operator==(const Image &other) const {
  return (width == other.width) && (height == other.height) &&
         (channels == other.channels) && (size == other.size) &&
         (pixels == other.pixels);
}

Image::Image(const Image &other)
    : channels{other.channels}, height{other.height}, width{other.width},
      size{other.size}, pixels{xt::xtensor<double, 3>(xt::zeros<double>(
                            {other.channels, other.height, other.width}))} {
  std::clog << "Copy Constructor\n";
  pixels = other.pixels;
}

Image::Image(Image &&other)
    : channels{other.channels}, height{other.height}, width{other.width},
      size{other.size}, pixels{std::move(other.pixels)} {
  std::clog << "Move Constructor\n";
}

// Image &Image::operator=(const Image &other) {
//   std::clog << "Copy Assignment Operator\n";
//   if (this != &other) {
//     channels = other.channels;
//     height = other.height;
//     width = other.width;
//     size = other.size;

//     pixels = xt::xtensor<double, 3>(
//         xt::zeros<double>({other.channels, other.height, other.width}));
//     pixels = other.pixels;
//   }
//   return *this;
// }

// Image &Image::operator=(Image &&other) {
//   std::clog << "Move Assignment Operator\n";
//   swap(other);
//   return *this;
// }

Image &Image::operator=(Image rhs) {
  std::clog << "Copy Assignment or Move Assignment Operator\n";
  swap(rhs);
  return *this;
}

void Image::swap(Image &other) {
  std::swap(channels, other.channels);
  std::swap(height, other.height);
  std::swap(width, other.width);
  std::swap(size, other.size);
  std::swap(pixels, other.pixels);
}

bool Image::save(std::string file_path) {
  /**
   * Save image as jpg or png file
   */
  auto file_extension = std::filesystem::path(file_path).extension();
  std::unique_ptr<unsigned char[]> out_data =
      std::make_unique<unsigned char[]>(width * height * channels);

  for (auto x = 0; x < width; ++x) {
    for (auto y = 0; y < height; ++y) {
      for (auto c = 0; c < channels; ++c) {
        int dst_idx = y * width * channels + x * channels + c;
        // Fill out_data with uint8 values range 0-255
        out_data[dst_idx] = std::roundf(pixels(c, y, x) * 255.0);
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
             file_extension == std::string(".PNG")) {
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

// void Image::set_pixel(int x, int y, int c, double val) {
//   if (x >= this->width || x < 0 || y >= this->height || y < 0 ||
//       c >= this->channels || c < 0) {
//     std::cerr << "set_pixel() error: Index out of bounds.\n";
//     std::exit(1);
//   }
//   this->data[c * this->width * this->height + y * this->width + x] = val;
// }

// double Image::get_pixel(int x, int y, int c) const {
//   if (x < 0)
//     x = 0;
//   if (x >= this->width)
//     x = this->width - 1;
//   if (y < 0)
//     y = 0;
//   if (y >= this->height)
//     y = this->height - 1;
//   return this->data[c * this->width * this->height + y * this->width + x];
// }

void Image::clamp() {
  // NOTE: np.where could be used here.
  for (auto x = 0; x < width; ++x) {
    for (auto y = 0; y < height; ++y) {
      for (auto c = 0; c < channels; ++c) {
        if (pixels(c, y, x) > 1.0) {
          pixels(c, y, x) = 1.0;
        }
        if (pixels(c, y, x) < 0.0) {
          pixels(c, y, x) = 0.0;
        }
      }
    }
  }
}

// map coordinate from 0-current_max range to 0-new_max range
double map_coordinate(double new_max, double current_max, double coord) {
  double a = new_max / current_max;
  double b = -0.5 + a * 0.5;
  return a * coord + b;
}

Image Image::resize(int new_w, int new_h, Interpolation method) const {
  Image resized(this->channels, new_h, new_w);
  double value = 0;
  for (int x = 0; x < new_w; x++) {
    for (int y = 0; y < new_h; y++) {
      for (int c = 0; c < resized.channels; c++) {
        double old_x = map_coordinate(this->width, new_w, x);
        double old_y = map_coordinate(this->height, new_h, y);
        if (method == Interpolation::BILINEAR) {
          value = bilinear_interpolate(*this, old_x, old_y, c);
        } else if (method == Interpolation::NEAREST) {
          ;
          value = nn_interpolate(*this, old_x, old_y, c);
        }
        // resized.set_pixel(x, y, c, value);
        resized.pixels(c, y, x) = value;
      }
    }
  }
  return resized;
}

double bilinear_interpolate(const Image &img, double x, double y, int c) {
  double p1, p2, p3, p4, q1, q2;
  double x_floor = std::floor(x), y_floor = std::floor(y);
  double x_ceil = x_floor + 1, y_ceil = y_floor + 1;
  p1 = img.pixels(c, y_floor, x_floor);
  p2 = img.pixels(c, y_floor, x_ceil);
  p3 = img.pixels(c, y_ceil, x_floor);
  p4 = img.pixels(c, y_ceil, x_ceil);
  q1 = (y_ceil - y) * p1 + (y - y_floor) * p3;
  q2 = (y_ceil - y) * p2 + (y - y_floor) * p4;
  return (x_ceil - x) * q1 + (x - x_floor) * q2;
}

double nn_interpolate(const Image &img, double x, double y, int c) {
  return img.pixels(c, std::round(y), std::round(x));
}

Image rgb_to_grayscale(const Image &img) {
  assert(img.channels == 3);
  Image gray(1, img.height, img.width);
  xt::xtensor<double, 2> red = xt::view(img.pixels, 0, xt::all(), xt::all());
  xt::xtensor<double, 2> green = xt::view(img.pixels, 1, xt::all(), xt::all());
  xt::xtensor<double, 2> blue = xt::view(img.pixels, 2, xt::all(), xt::all());
  xt::view(gray.pixels, 0, xt::all(), xt::all()) =
      0.299 * red + 0.587 * green + 0.114 * blue;
  return gray;
}

// Image grayscale_to_rgb(const Image &img) {
//   assert(img.channels == 1);
//   Image rgb(3, img.height, img.width);
//   xt::xtensor<double> red = xt::view(img.pixels, 0, xt::all(), xt::all());

//   return rgb;
// }

// separable 2D gaussian blur for 1 channel image
// Image gaussian_blur(const Image &img, double sigma) {
//   assert(img.channels == 1);

//   int size = std::ceil(6 * sigma);
//   if (size % 2 == 0)
//     size++;
//   int center = size / 2;
//   Image kernel(size, 1, 1);
//   double sum = 0;
//   for (int k = -size / 2; k <= size / 2; k++) {
//     double val = std::exp(-(k * k) / (2 * sigma * sigma));
//     kernel.set_pixel(center + k, 0, 0, val);
//     sum += val;
//   }
//   for (int k = 0; k < size; k++)
//     kernel.data[k] /= sum;

//   Image tmp(img.width, img.height, 1);
//   Image filtered(img.width, img.height, 1);

//   // convolve vertical
//   for (int x = 0; x < img.width; x++) {
//     for (int y = 0; y < img.height; y++) {
//       double sum = 0;
//       for (int k = 0; k < size; k++) {
//         int dy = -center + k;
//         sum += img.get_pixel(x, y + dy, 0) * kernel.data[k];
//       }
//       tmp.set_pixel(x, y, 0, sum);
//     }
//   }
//   // convolve horizontal
//   for (int x = 0; x < img.width; x++) {
//     for (int y = 0; y < img.height; y++) {
//       double sum = 0;
//       for (int k = 0; k < size; k++) {
//         int dx = -center + k;
//         sum += tmp.get_pixel(x + dx, y, 0) * kernel.data[k];
//       }
//       filtered.set_pixel(x, y, 0, sum);
//     }
//   }
//   return filtered;
// }

// void draw_point(Image &img, int x, int y, int size) {
//   for (int i = x - size / 2; i <= x + size / 2; i++) {
//     for (int j = y - size / 2; j <= y + size / 2; j++) {
//       if (i < 0 || i >= img.width)
//         continue;
//       if (j < 0 || j >= img.height)
//         continue;
//       if (std::abs(i - x) + std::abs(j - y) > size / 2)
//         continue;
//       if (img.channels == 3) {
//         img.set_pixel(i, j, 0, 1.f);
//         img.set_pixel(i, j, 1, 0.f);
//         img.set_pixel(i, j, 2, 0.f);
//       } else {
//         img.set_pixel(i, j, 0, 1.f);
//       }
//     }
//   }
// }

// void draw_line(Image &img, int x1, int y1, int x2, int y2) {
//   if (x2 < x1) {
//     std::swap(x1, x2);
//     std::swap(y1, y2);
//   }
//   int dx = x2 - x1, dy = y2 - y1;
//   for (int x = x1; x < x2; x++) {
//     int y = y1 + dy * (x - x1) / dx;
//     if (img.channels == 3) {
//       img.set_pixel(x, y, 0, 0.f);
//       img.set_pixel(x, y, 1, 1.f);
//       img.set_pixel(x, y, 2, 0.f);
//     } else {
//       img.set_pixel(x, y, 0, 1.f);
//     }
//   }
// }
