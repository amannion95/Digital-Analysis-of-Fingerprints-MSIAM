#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include "Picture.h"
using namespace cv;
using namespace std;








int main(int argc, char** argv )
{

  Mat image2 = imread( "../../data/blurred_finger.png" , IMREAD_COLOR );

  namedWindow("Display Image", WINDOW_NORMAL );

  Picture img(image2);

  cout << img.maximum_intensity() << " "  << img.minimum_intensity() << "\n";


  imshow("Display Image", image2);
    waitKey(0);
    return 0;

}
