#ifndef IMAGE_H
#define IMAGE_H

#include <string>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>

enum Interpolation { BILINEAR, NEAREST };

struct Image {
  explicit Image(const xt::xtensor<double, 3> &input_matrix);
  explicit Image(std::string file_path);
  Image(int w, int h, int c);
  Image();
  ~Image();

  Image(const Image &other);
  Image(Image &&other);

  Image &operator=(Image rhs);

  bool operator==(const Image &other) const;

  int width;
  int height;
  int channels;
  int size;
  xt::xtensor<double, 3> pixels;

  bool save(std::string file_path);
  void set_pixel(int x, int y, int c, double val);
  double get_pixel(int x, int y, int c) const;
  void clamp();
  Image resize(int new_w, int new_h, Interpolation method = BILINEAR) const;
  void _swap(Image &other);
};

double bilinear_interpolate(const Image &img, double x, double y, int c);
double nn_interpolate(const Image &img, double x, double y, int c);

Image rgb_to_grayscale(const Image &img);
Image grayscale_to_rgb(const Image &img);

Image gaussian_blur(const Image &img, double sigma);
double map_coordinate(double new_max, double current_max, double coord);

void draw_point(Image &img, int x, int y, int size = 3);
void draw_line(Image &img, int x1, int y1, int x2, int y2);

#endif
