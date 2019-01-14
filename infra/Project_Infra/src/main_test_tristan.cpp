#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include "Picture.h"
using namespace cv;
using namespace std;








int main(int argc, char** argv )
{

  Point c(100,150);
  Mat image(200,200,CV_8UC1,155);
  Mat image2 = imread( "/home/tristan/Bureau/Projet/Digital-Analysis-of-Fingerprints-MSIAM/infra/data/weak_finger.png" , 0 );
  //Mat image2 = imread( "/home/tristan/Bureau/Projet/Digital-Analysis-of-Fingerprints-MSIAM/infra/data/rtxy_small_finger.png" , 0 );
  namedWindow("Display Image", WINDOW_NORMAL );
  Picture img(image);
  Picture img2(image2);
  Picture img3("/home/tristan/Bureau/Projet/Digital-Analysis-of-Fingerprints-MSIAM/infra/data/moist_finger.png");
  Picture img4(image);

  //img2.print_picture();
  //img=img2.apply_threshold(0.1);
  //img4=img3.apply_threshold(0.1);
  img3=img2.apply_threshold(0.1);
  img3.apply_gaussian_blur(15);
  c=img3.pressure_center_gauss();



  img=img2.try_apply_isotrop(c,50,50);

  img2.print_picture();
  img.print_picture();
  //img.print_picture();
  //img.print_pression_center(101);
  //img4.print_pression_center(15);
  //img4.print_median_center(0.1);

  return 0;

}
