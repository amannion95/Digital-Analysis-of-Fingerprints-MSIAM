#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include "Picture.h"

using namespace cv;
using namespace std;

int main(int argc, char** argv )
{
  if(argc != 2){
    cerr << "usage: BoxImage <imagepath>" << endl;
    return EXIT_FAILURE;
  }

  Mat image = imread(argv[1], IMREAD_GRAYSCALE);

  /*
  unsigned int numboxes;
  cout << "Enter a number of boxes: " << endl;
  cin >> numboxes;
  */

  unsigned int input_intensity;
  cout << "Give a grayscale intensity for the boxes (between 0 and 255): " << endl;
  cin >> input_intensity;

  unsigned int x = image.rows;
  unsigned int y = image.cols;

  if(input_intensity < 0 || input_intensity > 255){
    cerr << "invalid grayscale intensity, closing" << endl;
    return EXIT_FAILURE;
  }

  else{
    for(unsigned int i=1; i<=4; i++){
      float x_float = x/4;
      float y_float = y/4;
      int x_size;
      int y_size;
      x_size = (int) y_float - 20;
      y_size = (int) x_float - 20;
      for(unsigned int j=(i-1)*x_size+10; j<=i*x_size; j++){
        for(unsigned int k=(i-1)*y_size+10; k<=i*y_size; k++){
          image.at<uchar>(j, k) = input_intensity;
        }
      }
    }
  }

  //saves image to data folder
  imwrite("../../data/boxesaddedimage.png", image);

  namedWindow("Display Image", WINDOW_NORMAL);
  imshow("Display Image", image);
  waitKey(0);
  return EXIT_SUCCESS;
}
