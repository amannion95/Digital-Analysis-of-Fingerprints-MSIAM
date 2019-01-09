#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include "Picture.h"
using namespace cv;
using namespace std;








int main(int argc, char** argv )
{
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
    return 0;

}
