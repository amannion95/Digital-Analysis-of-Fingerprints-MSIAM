#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include "Picture.h"
using namespace cv;
using namespace std;








int main(int argc, char** argv )
{

  Mat image2 = imread( "../../data/blurred_finger.png" , IMREAD_COLOR );



    Mat image(100,100,CV_8UC1,21);

    namedWindow("Display Image", WINDOW_NORMAL );



    for(int j=0;j<=100  ;j++){
      for(int i=100;i<=199;i++){
        image2.at<uchar>(j,i)=0;

}
}
  imshow("Display Image", image2);
    waitKey(0);
    return 0;

}
