#include <cstdlib>
#include <iostream>
#include <vector>
#include <string.h>
#include <opencv2/opencv.hpp>

#include "Picture.h"
#include "Useful_functions.h"

int main(){
  Mat image = imread("../../data/clean_finger.png");
  Picture img(image);

  Point p(90, 220);

  Picture 4zones = img.swirl_zonalmorph(p, 6.5, 80, 4);
  4zones.print_picture("Swirl function with zonal morphology: 4x4 grid");

  Picture 5zones = img.swirl_zonalmorph(p, 6.5, 80, 5);
  5zones.print_picture("Swirl function with zonal morphology: 5x5 grid");

  Picture 10zones = img.print_picture(p, 6.5, 80, 10);
  10zones.print_picture("Swirl function with zonal morphology: 10x10 grid");

  return 0;
}
