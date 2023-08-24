#include <iostream> 
#include <string>

#include "image.hpp"
#include "sift.hpp"

int main(int argc, char *argv[])
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    if (argc != 2) {
        std::cerr << "Usage: ./find_keypoints input.jpg (or .png)\n";
        return 0;
    }
    // auto tmp =Eigen::MatrixXd(2,3);
    // tmp << 1,2,3,4,5,6;
    Image img(argv[1]);
    img =  img.channels == 1 ? img : rgb_to_grayscale(img);

    std::cout << img.width << img.height << img.size << std::endl;

    std::vector<sift::Keypoint> kps = sift::find_keypoints_and_descriptors(img);
    Image result = sift::draw_keypoints(img, kps);
    result.save("find_keypoints_result.png");

    std::cout << "Found " << kps.size() << " keypoints. Output image is saved as a png file.\n";
    return 0;
}
