#include <gtest/gtest.h>
#include <image.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <xtensor/xadapt.hpp>

#include <limits.h>
#include <stdio.h>
#include <unistd.h>

TEST(TestFindCurrentPath, TestFindCurrentPath) {
  char cwd[PATH_MAX];
  printf("================================================================\n");
  if (getcwd(cwd, sizeof(cwd)) != NULL) {
    printf("Current working dir: %s\n", cwd);
  } else {
    perror("getcwd() error");
  }
  printf("================================================================\n");
}

Image get_ones_image(int channel, int height, int width) {
  auto img = Image(channel, height, width);
  img.pixels += 1;
  return img;
}

TEST(TestImage, map_coordinate) {
  // 0-99: 100
  // 0-24: 25
  double cur_max = 99.0;
  double new_max = 24.0;
  double coord = 99.0;
  double new_coord = map_coordinate(new_max, cur_max, coord);
  assert(int(new_coord) == 24);
  coord = 0.0;
  new_coord = map_coordinate(new_max, cur_max, coord);
  assert(int(new_coord) == 0);

  cur_max = 24.0;
  new_max = 99.0;
  coord = 24.0;
  new_coord = map_coordinate(new_max, cur_max, coord);
  assert(int(new_coord) == 99);
  coord = 0.0;
  new_coord = map_coordinate(new_max, cur_max, coord);
  assert(int(new_coord) == 0);
}

TEST(TestImage, TestWithArray) {
  Image test_img(2, 3, 4);
  int counter = 0;
  for (int x = 0; x < test_img.width; x++) {
    for (int y = 0; y < test_img.height; y++) {
      for (int c = 0; c < test_img.channels; c++) {
        EXPECT_EQ(test_img.pixels(c, y, x), 0.0);
        test_img.pixels(c, y, x) = counter;
        ++counter;
      }
    }
  }
  std::clog << "test_img:\n" << test_img.pixels << '\n';

  std::clog << "Test Copy Constructor.\n";
  Image img_cp_ctor{test_img};
  std::clog << "img_cp_ctor:\n" << img_cp_ctor.pixels << '\n';
  EXPECT_EQ(img_cp_ctor, test_img);

  std::clog << "Test Copy Assignment Operator.\n";
  Image img_cp_asgn_otor;
  img_cp_asgn_otor = test_img;
  std::clog << "img_cp_asgn_otor\n" << img_cp_asgn_otor.pixels << '\n';
  EXPECT_EQ(img_cp_asgn_otor, test_img);
  EXPECT_EQ(img_cp_ctor, test_img);

  std::clog << "Test Move Assignment Operator.\n";
  Image img_mv_asgn_otor;
  img_mv_asgn_otor = std::move(test_img);
  // test_img becomes unspecified after it's moved from.
  EXPECT_EQ(img_mv_asgn_otor, img_cp_ctor);

  std::clog << "Test Move Constructor.\n";
  Image img_mv_ctor{get_ones_image(2, 3, 4)};
  std::clog << "img_mv_ctor:\n" << img_mv_ctor.pixels << '\n';

  std::vector<Image> vec;
  vec.push_back(get_ones_image(2, 3, 4));
  std::clog << "vec[0]:\n" << vec[0].pixels << '\n';
  Image tmp_vec = get_ones_image(2, 3, 4);
  EXPECT_EQ(tmp_vec, vec[0]);
}

TEST(TestImage, JPGImageReadTest) {
  std::string jpg_path = "../../../imgs/book_in_scene.jpg";
  Image test_img(jpg_path);
  EXPECT_EQ(test_img.channels, 3);
  EXPECT_EQ(test_img.height, 540);
  EXPECT_EQ(test_img.width, 720);
}

TEST(TestImage, PNGImageReadTest) {
  std::string png_path = "../../../imgs/book.png";
  Image test_img(png_path);
  EXPECT_EQ(test_img.channels, 3);
  EXPECT_EQ(test_img.height, 556);
  EXPECT_EQ(test_img.width, 394);
}

TEST(TestImage, JPGImageWriteTest) {
  std::string img_path = "../../../imgs/book_in_scene.jpg";
  std::string img_output_path =
      "../../../tests/tmp/image/book_in_scene_gray.jpg";

  Image test_img(img_path);
  Image gray = rgb_to_grayscale(test_img);
  gray.save(img_output_path);
}

TEST(TestImage, PNGImageWriteTest) {
  std::string img_path = "../../../imgs/book.png";
  std::string img_output_path = "../../../tests/tmp/image/book_gray.png";

  Image test_img(img_path);
  Image gray = rgb_to_grayscale(test_img);
  gray.save(img_output_path);
}

TEST(TestImage, PNGResizeImage) {
  std::string img_path = "../../../imgs/book.png";
  std::string img_output_path = "../../../tests/tmp/image/book_resize.png";

  Image test_img(img_path);

  Image half_img = test_img.resize(int(test_img.width / 5),
                                   int(test_img.height / 5), BILINEAR);

  half_img.save(img_output_path);
}

TEST(TestImage, JPGResizeImage) {
  std::string img_path = "../../../imgs/book_in_scene.jpg";
  std::string img_output_path =
      "../../../tests/tmp/image/book_in_scene_resize.jpg";

  Image test_img(img_path);

  Image half_img = test_img.resize(int(test_img.width / 5),
                                   int(test_img.height / 5), BILINEAR);

  half_img.save(img_output_path);
}

TEST(TestImage, PNGGaussianBlurImage) {
  std::string img_path = "../../../imgs/book.png";
  std::string img_output_path =
      "../../../tests/tmp/image/book_gaussianblur.png";

  Image test_img(img_path);
  Image gray = rgb_to_grayscale(test_img);

  double sigma = 2.0;
  Image gauss = gaussian_blur(gray, sigma);
  gauss.save(img_output_path);
}

TEST(TestImage, JPGGaussianBlurImage) {
  std::string img_path = "../../../imgs/book_in_scene.jpg";
  std::string img_output_path =
      "../../../tests/tmp/image/book_in_scene_gaussianblur.jpg";

  Image test_img(img_path);
  Image gray = rgb_to_grayscale(test_img);

  double sigma = 2.0;
  Image gauss = gaussian_blur(gray, sigma);
  gauss.save(img_output_path);
}

TEST(TestImage, PNGGaussianBlurAndResizeImage) {
  std::string img_path = "../../../imgs/book.png";
  std::string img_output_path =
      "../../../tests/tmp/image/book_gaussianblur_resize.png";

  // TODO FIND OUT THE BETTER WAY TO DO THIS
  // red = gaussian_blur(test_img[0, :, :])
  Image test_img(img_path);
  Image gray = rgb_to_grayscale(test_img);

  double sigma = 2.5;
  Image gauss = gaussian_blur(gray, sigma);
  Image half_img =
      gauss.resize(int(gauss.width / 5), int(gauss.height / 5), BILINEAR);

  Image rgb = grayscale_to_rgb(half_img);
  assert(rgb.channels == 3);

  rgb.save(img_output_path);
}

TEST(TestImage, JPGGaussianBlurAndResizeImage) {
  std::string img_path = "../../../imgs/book_in_scene.jpg";
  std::string img_output_path =
      "../../../tests/tmp/image/book_in_scene_gaussianblur_resize.jpg";

  // TODO FIND OUT THE BETTER WAY TO DO THIS
  // red = gaussian_blur(test_img[0, :, :])
  Image test_img(img_path);
  Image gray = rgb_to_grayscale(test_img);

  double sigma = 2.5;
  Image gauss = gaussian_blur(gray, sigma);
  Image half_img =
      gauss.resize(int(gauss.width / 5), int(gauss.height / 5), BILINEAR);

  Image rgb = grayscale_to_rgb(half_img);
  assert(rgb.channels == 3);

  rgb.save(img_output_path);
}
