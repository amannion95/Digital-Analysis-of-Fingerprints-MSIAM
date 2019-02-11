#include <cstdlib>
#include <iostream>
#include <opencv2/opencv.hpp>

#include "Picture.h"
#include "Useful_functions.h"

int main(){
  Mat image = imread("../../data/clean_finger.png", IMREAD_GRAYSCALE);
  Picture img(image);

  Picture img1 = img.swirl(Point(90,220), 6.5, 80);
  img1.print_picture("Locally warped image (geometry only)");

  Picture img2 = img1.local_erosion(Point(90,210), 60, 90);
  Picture img3 = img2.local_dilation(Point(220,260), 100, 100);
  img3.print_picture("Locally warped image with morphological effects");

  return 0;
}
