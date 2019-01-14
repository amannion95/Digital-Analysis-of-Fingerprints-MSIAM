#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include "Picture.h"

using namespace cv;
using namespace std;

int main(int argc, char** argv )
{
  Mat image =  imread("../../data/clean_finger.png", IMREAD_GRAYSCALE);

  Picture img(image);
  Point p = img.pressure_center_gauss();

  vector<Point> v = img.ellipse_nbh(p, 30, 50);
  img.show_nbh(v);

  Picture pimg = img.log_transform_isotropic(p, 30, 50, 0.1);

  pimg.print_picture();
  /*
  Mat image(200,200,CV_8UC1,125);
  Mat image2 = imread( "../../data/blurred_finger.png" , IMREAD_COLOR );

  namedWindow("Display Image", WINDOW_NORMAL );

  Picture img(image);
  Picture img2(image2);

  cout << img.maximum_intensity() << " "  << img.minimum_intensity() << "\n";
  cout << img2.maximum_intensity() << " "  << img2.minimum_intensity() << "\n";

  for(int i=0;i<50;i++){
    for(int j=0;j<50;j++){
      image.at<uchar>(j,i)=200;
    }
  }
  imshow("Display Image", image);
  waitKey(0);
  Picture img_square(image);
  cout << img_square.maximum_intensity() << " "  << img_square.minimum_intensity() << "\n";
  img_square.rescale_color();
  img_square.print_picture();
    waitKey(0);

    /*Mat image2 = imread( "../../data/squares.png" , IMREAD_GRAYSCALE );
    //cvtColor( image, image2, CV_BGR2GRAY );
    namedWindow("Display Image", WINDOW_NORMAL );
<<<<<<< HEAD
    Picture pic(image2);
    Mat img (300,200,CV_8UC1);
    Picture img_pic(img);
=======


    Picture pic(image2);

    Mat img (300,200,CV_8UC1);
    Picture img_pic(img);

>>>>>>> f87c428abb79fb9386c127916096ca74321b55e1
    for(int j=0;j<img_pic.get_y_len()-50 ;j++){
      for(int i=0;i<img_pic.get_x_len()-50;i++){
        img_pic.set_intensity(j,i,1);
    }
  }
  //pic.print_picture();
  //pic.symmetry_wrt_y().print_picture();
  //pic.symmetry_wrt_x().print_picture();
  //img_pic.print_picture();
  //img_pic.symmetry_wrt_y().print_picture();
  //img_pic.diagonal_symmetry_top_to_bottom().print_picture();
  //img_pic.diagonal_symmetry_bottom_to_top().print_picture();
  */
  return 0;

}
