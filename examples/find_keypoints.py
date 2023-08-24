import lion_sift_cpp
import numpy as np
from dataclasses import dataclass
from PIL import Image
import sys

""" The returned value of lion_sift_cpp.computeKeypointsAndDescriptors()
is a data structure defined below:

struct Keypoint {
  // discrete coordinates
  int i;
  int j;
  int octave;
  int scale; // index of gaussian image inside the octave

  // continuous coordinates (interpolated)
  double x;
  double y;
  double sigma;
  double extremum_val; // value of interpolated DoG extremum

  std::array<uint8_t, 128> descriptor;
};
"""


if __name__ == "__main__":
    # Open a PNG file
    input_path = sys.argv[1]
    print(input_path)
    pil_img = Image.open(input_path)

    # Convert the image to a NumPy array
    # pixels: (height, width, channel)
    pixels = np.asarray(pil_img)[:, :, 0:1]

    print(f"Input image's shape: {pixels.shape}.") 

    # `pixels` must be (channel, height, width) to be as input.
    pixels = np.moveaxis(pixels, -1, 0)
    print(f"Convert shape to {pixels.shape}.")
    keypoints_and_descriptors = lion_sift_cpp.computeKeypointsAndDescriptors(pixels)

    print(f"Found {len(keypoints_and_descriptors)} keypoints.")
    for kp in keypoints_and_descriptors:
        print(f"The first keypoint:")
        print(f"kp.i: {kp.i}, kp.j: {kp.j}, kp.x: {kp.x}, kp.y: {kp.y}")
        descriptor = np.array(kp.descriptor)
        print(f"kp.descriptor is a 128-vector whose shape is {descriptor.shape}.")
        break