#include <cstdlib>
#include <iostream>
#include <opencv2/opencv.hpp>

#include "Picture.h"
#include "Useful_functions.h"

int main(){
  Mat image = imread("../../data/clean_finger.png", IMREAD_GRAYSCALE);
  Picture img(image);

  Point p(128, 144);

  Picture img1 = img.cast_rotation_polar(p, 45.0);
  Picture err1 = img1 - img;
  img1.print_picture("Rotation using casted intensities");
  err1.print_picture("Error image for casted intensities");

  Picture img2 = img.nn_rotation_polar(p, 45.0);
  Picture err2 = img2 - img;
  img2.print_picture("Rotation using nearest-neighbour interpolation");
  err2.print_picture("Error image for nearest-neighbour interpolation");

  Picture img3 = img.bilinear_rotation_polar(p, 45.0);
  Picture err3 = img3 - img;
  img3.print_picture("Rotation using bilinear interpolation");
  err3.print_picture("Error image for bilinear interpolation");

  return 0;
}
