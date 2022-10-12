from build import pybind11_wrapper
import numpy as np
import pickle
from typing import List, Tuple, Union
import rasterio
import rpcm
from dataclasses import dataclass
from rasterio.windows import Window

@dataclass
class Keypoint:
    i: int # x coord
    j: int # y coord
    octave: int
    scale: int #index of gaussian image inside the octave

    # continuous coordinates (interpolated)
    x: float
    y: float
    sigma: float
    extremum_value: float # value of interpolated DoG extremum
    
def transform_keypoint_coordinate(
    list_keypoint, offset: Tuple[int, int]
):

    transformed_keypoint = [
        Keypoint(
            i=kp.i + offset[0],
            j=kp.j + offset[1],
            octave=kp.octave,
            scale=kp.scale,
            x=kp.x + offset[0],
            y=kp.y + offset[1],
            sigma=kp.sigma,
            extremum_value=kp.extremum_val,
        )
        for kp in list_keypoint
    ]
    descriptor = [kp.descriptor for kp in list_keypoint]
    return transformed_keypoint, np.array(descriptor)

if __name__ == "__main__":
    img1 = "/home/lai/LiveEO/09_photogrammetry/tests/test_pleiades/IMG_PHR1A_P_201807271102180_SEN_3385144101-1_R1C1.JP2"
    rpc1 = rpcm.rpc_from_rpc_file(
        "/home/lai/LiveEO/09_photogrammetry/tests/test_pleiades/RPC_PHR1A_P_201807271102180_SEN_3385144101-1.XML"
    )
    img2 = "/home/lai/LiveEO/09_photogrammetry/tests/test_pleiades/IMG_PHR1A_P_201807271101330_SEN_3385144101-1_R1C1.JP2"
    rpc2 = rpcm.rpc_from_rpc_file(
        "/home/lai/LiveEO/09_photogrammetry/tests/test_pleiades/RPC_PHR1A_P_201807271101330_SEN_3385144101-1.XML"
    )
    x1 = 7195
    y1 = 7432
    w1 = 500
    h1 = 500
    with rasterio.open(img1, "r") as src:
        crop1 = src.read(1, window=Window(x1, y1, w1, h1))

    x2 = 7513
    y2 = 7580
    w2 = 565
    h2 = 665
    with rasterio.open(img2, "r") as src:
        crop2 = src.read(1, window=Window(x2, y2, w2, h2))
    print(crop1.shape, crop2.shape)
    keypoints_and_descriptors1 = pybind11_wrapper.computeKeypointsAndDescriptors(crop1)
    keypoints_and_descriptors2 = pybind11_wrapper.computeKeypointsAndDescriptors(crop2)
    for kp in keypoints_and_descriptors1:
        print(f"kp.i: {kp.i}, kp.j: {kp.j}, kp.x: {kp.x}, kp.y: {kp.y}")
        break
    keypoints1, descriptors1 = transform_keypoint_coordinate(keypoints_and_descriptors1, (x1, y1))
    keypoints2, descriptors2 = transform_keypoint_coordinate(keypoints_and_descriptors2, (x2, y2))

    # for kp in keypoints1:
    #     print(f"kp.i: {kp.i}, kp.j: {kp.j}, kp.x: {kp.x}, kp.y: {kp.y}")
    print(f"Image 1: {len(keypoints1)}")
    print(f"Image 2: {len(keypoints2)}")

    with open("/home/lai/LiveEO/09_photogrammetry/cpp_sift_keypoints_descriptors_left.pkl", "wb") as output_pickle:
        temp = {}
        pickle_keypoints = []
        for kp in keypoints1:
            pickle_keypoints.append(
                (
                    kp.x, kp.y
                )
            )
        temp.update(keypoints=pickle_keypoints)
        temp.update(descriptor=descriptors1)
        pickle.dump(temp, output_pickle, pickle.HIGHEST_PROTOCOL)
    with open("/home/lai/LiveEO/09_photogrammetry/cpp_sift_keypoints_descriptors_right.pkl", "wb") as output_pickle:
        temp = {}
        pickle_keypoints = []
        for kp in keypoints2:
            pickle_keypoints.append(
                (
                    kp.x, kp.y
                )
            )
        temp.update(keypoints=pickle_keypoints)
        temp.update(descriptor=descriptors2)
        pickle.dump(temp, output_pickle, pickle.HIGHEST_PROTOCOL)
        
    with open("/home/lai/LiveEO/09_photogrammetry/cpp_sift_keypoints_descriptors_left.pkl", "rb") as input_pickle:
        input_pickle_obj = pickle.load(input_pickle)
        print(len(input_pickle_obj["keypoints"]))
        print(input_pickle_obj["descriptor"].shape)

        for index, keypoint in enumerate(input_pickle_obj["keypoints"]):
            #print(keypoints1[index])
            assert keypoint[0] < x1+w1
            assert keypoint[0] > x1
            assert keypoint[1] < y1+h1
            assert keypoint[1] > y1
        #print(input_pickle_obj["descriptor"][0,:])
    with open("/home/lai/LiveEO/09_photogrammetry/cpp_sift_keypoints_descriptors_right.pkl", "rb") as input_pickle:
        input_pickle_obj = pickle.load(input_pickle)
        print(len(input_pickle_obj["keypoints"]))
        print(input_pickle_obj["descriptor"].shape)
        for keypoint in input_pickle_obj["keypoints"]:
            assert keypoint[0] < x2+w2
            assert keypoint[0] > x2
            assert keypoint[1] < y2+h2
            assert keypoint[1] > y2
