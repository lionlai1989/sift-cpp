#include <gtest/gtest.h>
#include <image.hpp>
#include <iostream>
#include <sift.hpp>
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

TEST(TestSIFT, PNGfind_keypoints_and_descriptors) {
  std::string img_path = "../../../imgs/book.png";
  std::string img_output_path = "../../../tests/tmp/sift/book_keypoints.png";

  Image test_img(img_path);
  std::vector<sift::Keypoint> kps =
      sift::find_keypoints_and_descriptors(test_img);
  Image result = sift::draw_keypoints(test_img, kps);
  result.save(img_output_path);
  std::cout << "Found " << kps.size() << " keypoints.";
  std::cout << "Output image is saved as " << img_output_path << '\n';
}

TEST(TestSIFT, JPGfind_keypoints_and_descriptors) {
  std::string img_path = "../../../imgs/book_in_scene.jpg";
  std::string img_output_path =
      "../../../tests/tmp/sift/book_in_scene_keypoints.jpg";

  Image test_img(img_path);
  std::vector<sift::Keypoint> kps =
      sift::find_keypoints_and_descriptors(test_img);
  Image result = sift::draw_keypoints(test_img, kps);
  result.save(img_output_path);
  std::cout << "Found " << kps.size() << " keypoints.";
  std::cout << "Output image is saved as " << img_output_path << '\n';
}

TEST(TestSIFT, find_keypoint_matches) {
  std::string img1_path = "../../../imgs/book.png";
  std::string img2_path = "../../../imgs/book_in_scene.jpg";
  std::string img_output_path =
      "../../../tests/tmp/sift/book_matching_kpts.png";

  Image img1(img1_path);
  Image img2(img2_path);

  std::vector<sift::Keypoint> kps_a =
      sift::find_keypoints_and_descriptors(img1);
  std::vector<sift::Keypoint> kps_b =
      sift::find_keypoints_and_descriptors(img2);
  std::vector<std::pair<int, int>> matches =
      sift::find_keypoint_matches(kps_a, kps_b);
  Image result = sift::draw_matches(img1, img2, kps_a, kps_b, matches);
  result.save(img_output_path);
  std::cout << "Found " << matches.size() << " matching points.";
  std::cout << "Output image is saved as " << img_output_path << '\n';
}

TEST(TestSIFT, RotatedBook) {
  std::string img1_path = "../../../imgs/book_rotated.jpg";
  std::string img2_path = "../../../imgs/book_in_scene.jpg";
  std::string img_output_path =
      "../../../tests/tmp/sift/book_matching_kpts.png";

  Image img1(img1_path);
  Image img2(img2_path);

  std::vector<sift::Keypoint> kps_a =
      sift::find_keypoints_and_descriptors(img1);
  std::vector<sift::Keypoint> kps_b =
      sift::find_keypoints_and_descriptors(img2);
  std::vector<std::pair<int, int>> matches =
      sift::find_keypoint_matches(kps_a, kps_b);
  Image result = sift::draw_matches(img1, img2, kps_a, kps_b, matches);
  result.save(img_output_path);
  std::cout << "Found " << matches.size() << " matching points.";
  std::cout << "Output image is saved as " << img_output_path << '\n';
}

TEST(TestSIFT, Box) {
  std::string img1_path = "../../../imgs/box.png";
  std::string img2_path = "../../../imgs/box_in_scene.png";
  std::string img_output_path = "../../../tests/tmp/sift/box_matching_kpts.png";

  Image img1(img1_path);
  Image img2(img2_path);

  std::vector<sift::Keypoint> kps_a =
      sift::find_keypoints_and_descriptors(img1);
  std::vector<sift::Keypoint> kps_b =
      sift::find_keypoints_and_descriptors(img2);
  std::vector<std::pair<int, int>> matches =
      sift::find_keypoint_matches(kps_a, kps_b);
  Image result = sift::draw_matches(img1, img2, kps_a, kps_b, matches);
  result.save(img_output_path);
  std::cout << "Found " << matches.size() << " matching points.";
  std::cout << "Output image is saved as " << img_output_path << '\n';
}

TEST(TestSIFT, Mountain) {
  std::string img1_path = "../../../imgs/mountain1.png";
  std::string img2_path = "../../../imgs/mountain2.png";
  std::string img_output_path =
      "../../../tests/tmp/sift/mountain_matching_kpts.png";

  Image img1(img1_path);
  Image img2(img2_path);

  std::vector<sift::Keypoint> kps_a =
      sift::find_keypoints_and_descriptors(img1);
  std::vector<sift::Keypoint> kps_b =
      sift::find_keypoints_and_descriptors(img2);
  std::vector<std::pair<int, int>> matches =
      sift::find_keypoint_matches(kps_a, kps_b);
  Image result = sift::draw_matches(img1, img2, kps_a, kps_b, matches);
  result.save(img_output_path);
  std::cout << "Found " << matches.size() << " matching points.";
  std::cout << "Output image is saved as " << img_output_path << '\n';
}

TEST(TestSIFT, Brush) {
  std::string img1_path = "../../../imgs/brush1.png";
  std::string img2_path = "../../../imgs/brush2.png";
  std::string img_output_path =
      "../../../tests/tmp/sift/brush_matching_kpts.png";

  Image img1(img1_path);
  Image img2(img2_path);

  std::vector<sift::Keypoint> kps_a =
      sift::find_keypoints_and_descriptors(img1);
  std::vector<sift::Keypoint> kps_b =
      sift::find_keypoints_and_descriptors(img2);
  std::vector<std::pair<int, int>> matches =
      sift::find_keypoint_matches(kps_a, kps_b);
  Image result = sift::draw_matches(img1, img2, kps_a, kps_b, matches);
  result.save(img_output_path);
  std::cout << "Found " << matches.size() << " matching points.";
  std::cout << "Output image is saved as " << img_output_path << '\n';
}
