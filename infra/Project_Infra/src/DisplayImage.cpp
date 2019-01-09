#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <cstdlib>

using namespace std;
using namespace cv;

int main(int argc, const char* argv[]){
  if(argc != 2){
    cerr << "usage: DisplayImage.out <ImagePath>" << endl;
    return EXIT_FAILURE;
  }

  //matrix data type
  Mat image;
  image = imread(argv[1], 1);

  if(!image.data){
    cerr << "no image data" << endl;
    return EXIT_FAILURE;
  }

  //create a window
  namedWindow("DisplayImage", WINDOW_AUTOSIZE);
  //load image in window
  imshow("DisplayImage", image);

  //without this the window will close straightaway
  waitKey(0);

  return EXIT_SUCCESS;
}
