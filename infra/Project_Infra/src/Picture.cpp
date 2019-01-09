#include "Picture.h"
#include "Usefull_functions.h"

using namespace cv;




Picture::Picture(const std::string& filename){
  picture=imread(filename);
  x_length=(picture.size()).height;
  y_length=(picture.size()).width;
}

Picture::Picture(unsigned int x_length=0,unsigned int y_length=0){
  Mat image(y_length,x_length,CV_8UC1);
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
  Picture symmetry(x_length,y_length);
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      symmetry.set_intensity(y_length-1-j,x_length-1-i,this->get_intensity(j,i));
    }
  }
  return symmetry;
}

void Picture::operator=(Picture Pic){
  picture=(Pic.picture).clone();
  x_length=Pic.x_length;
  y_length=Pic.y_length;
}
