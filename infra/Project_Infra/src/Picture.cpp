#include "Picture.h"
#include "Usefull_functions.h"

using namespace cv;
using namespace std;



Picture::Picture(const std::string& filename){
  picture=cv::imread(filename);
  x_length=(picture.size()).height;
  y_length=(picture.size()).width;
}


Picture::Picture(unsigned int x_length,unsigned int y_length){
  Mat image(y_length,x_length,CV_8UC1 );
}

Picture::Picture(const cv::Mat& pic){
  picture=pic.clone();
  x_length=(picture.size()).width;
  y_length=(picture.size()).height;
}

float Picture::get_intensity(unsigned int i, unsigned int j)const{
  return iitof((int)picture.at<uchar>(i,j));
}

void Picture::set_intensity(unsigned int i, unsigned int j,float intensity){
  if ((intensity<0)||(intensity>1)){
    std::cerr<<"Wrong intensity value, she must belong to [0,1]"<<std::endl;
  }
  picture.at<uchar>(i,j)=iftoi(intensity);
}

void Picture::print_picture()const{
  namedWindow("Display Image", WINDOW_NORMAL );
  imshow("Display Image", picture);
  waitKey(0);
}

float Picture::maximum_intensity(){
  float max_intensity=0.;
    for (int i = 0 ; i < x_length ; i++ ){
      for (int j = 0 ; j < y_length ; j++ ){
        max_intensity=max(get_intensity(i,j),max_intensity);
      }
    }
   return(max_intensity);
 }

float Picture::minimum_intensity(){
  float min_intensity=1.;
    for (int i = 0 ; i < x_length ; i++ ){
      for (int j = 0 ; j < y_length ; j++ ){
        min_intensity=min(get_intensity(i,j),min_intensity);
      }
    }
  return(min_intensity);
}
