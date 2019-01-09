#include "Picture.h"
#include "Usefull_functions.h"

using namespace cv;




Picture::Picture(const std::string& filename){
  picture=imread(filename);
  x_length=(picture.size()).height;
  y_length=(picture.size()).width;
}

Picture::Picture(unsigned int x_length,unsigned int y_length){
  Mat image(y_length,x_length,CV_8UC1);

  //check if it is y_len then x_len or x_len then y_len

  picture=image.clone();
  this->x_length=x_length;
  this->y_length=y_length;
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


Picture Picture::symmetry_transform()const{
  Picture tmp(x_length,y_length);
  Picture symmetry(x_length,y_length);
  if (x_length%2==0){
    for(int i=0;i<x_length/2;i++){
      for(int j=0;j<y_length;j++){
        tmp.set_intensity(i,j,)
    }
  }
  else {

  }


  return symmetry;
}
