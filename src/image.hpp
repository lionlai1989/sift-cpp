#ifndef IMAGE_H
#define IMAGE_H
#include <string>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
enum Interpolation { BILINEAR, NEAREST };

struct Image {
  Image();
  explicit Image(const xt::xtensor<double, 3> &input_matrix);
  explicit Image(std::string file_path);
  Image(int c, int h, int w);

  Image(const Image &other);
  Image(Image &&other);

  //   Image &operator=(const Image &other);
  //   Image &operator=(Image &&other);
  Image &operator=(Image rhs);
  ~Image();

  bool operator==(const Image &other) const;

  int channels;
  int height;
  int width;
  int size;
  xt::xtensor<double, 3> pixels;

  bool save(std::string file_path);
  void swap(Image &other);

  void clamp();
  Image resize(int new_w, int new_h, Interpolation method = BILINEAR) const;
};

double bilinear_interpolate(const Image &img, double x, double y, int c);
double nn_interpolate(const Image &img, double x, double y, int c);

Image rgb_to_grayscale(const Image &img);
// Image grayscale_to_rgb(const Image &img);

// Image gaussian_blur(const Image &img, double sigma);

// void draw_point(Image &img, int x, int y, int size = 3);
// void draw_line(Image &img, int x1, int y1, int x2, int y2);

#endif
