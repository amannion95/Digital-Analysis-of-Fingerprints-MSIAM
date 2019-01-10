#include "Picture.h"
#include "Usefull_functions.h"
#include<opencv2/opencv.hpp>
#include <math.h>
using namespace cv;




Picture::Picture(const std::string& filename){
  picture=imread(filename, IMREAD_GRAYSCALE);
  x_length=(picture.size()).height;
  y_length=(picture.size()).width;
}

Picture::Picture(unsigned int x_length,unsigned int y_length){
  Mat image(y_length,x_length,CV_8UC1);
  picture=image.clone();
  this->x_length=x_length;
  this->y_length=y_length;
}

Picture::Picture(const cv::Mat& pic){
  picture=pic.clone();
  x_length=(picture.size[1]);
  y_length=(picture.size[0]);
}

Picture::Picture(){
  x_length=0;
  y_length=0;
  Mat image(0,0,CV_8UC1);
  picture=image;
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
  std::cout<<"X_len : "<<x_length<<" Y_len : "<<y_length<<std::endl;
  namedWindow("Display Image", WINDOW_NORMAL );
  imshow("Display Image", picture);
  waitKey(0);
}


Picture Picture::symmetry_wrt_y()const{
  Picture symmetry;
  symmetry=clone();
  std::cout<<symmetry.x_length<<std::endl;
  for(int j=0;j<y_length;j++){
    for(int i=0;i<x_length;i++){
      symmetry.set_intensity(j,x_length-1-i,get_intensity(j,i));
    }
  }

  return symmetry;
}


Picture Picture::symmetry_wrt_x()const{
  Picture symmetry;
  symmetry=clone();
  for(int j=0;j<y_length;j++){
    for(int i=0;i<x_length;i++){
      symmetry.set_intensity(y_length-1-j,i,get_intensity(j,i));
    }
  }

  return symmetry;
}
Picture Picture::diagonal_symmetry_top_to_bottom()const{
  Picture sym(picture.t());
  return sym;
}



Picture Picture::diagonal_symmetry_bottom_to_top()const{
  Picture sym(picture.t());
  return sym.symmetry_wrt_y().symmetry_wrt_x();
}

void Picture::operator=(Picture Pic){
  picture=(Pic.picture).clone();
  x_length=Pic.x_length;
  y_length=Pic.y_length;
}

float Picture::maximum_intensity()const{
  float max_intensity=0.;
    for (int i = 0 ; i < x_length ; i++ ){
      for (int j = 0 ; j < y_length ; j++ ){
        max_intensity=max(get_intensity(i,j),max_intensity);
      }
    }
   return(max_intensity);
}

float Picture::minimum_intensity()const{
  float min_intensity=1.;
    for (int i = 0 ; i < x_length ; i++ ){
      for (int j = 0 ; j < y_length ; j++ ){
        min_intensity=min(get_intensity(i,j),min_intensity);
      }
    }
  return(min_intensity);
}

Picture Picture::clone()const{
  Picture clone;
  clone.picture=picture.clone();
  clone.x_length=x_length;
  clone.y_length=y_length;
  return clone;
}


unsigned int Picture::get_x_len(){
  return x_length;
}
unsigned int Picture::get_y_len(){
  return y_length;
}

void Picture::rescale_color(){
  float min_intensity = minimum_intensity();
  float max_intensity = maximum_intensity();
  float size = max_intensity - min_intensity;
  assert ( size != 0 ) ;
  for (int i = 0 ; i < x_length ; i++ ){
    for (int j = 0 ; j < y_length ; j++ ){
      set_intensity(i,j,(get_intensity(i,j)-min_intensity)/size);
    }
  }
}

float** Picture::get_matrix(){
  int row=x_length;
  int col=y_length;
  float **matrix = new float*[row];
  for ( int i = 0 ; i < row ; i ++ ) {
    matrix[i] = new float[col];
  }
  for(int i = 0 ; i < row; i++ ){
    for(int j = 0 ; j < col ; j++ ){
      matrix[i][j] = iitof(picture.data[i*row+j]);
    }
  }
  return matrix;
}
float** Picture::apply_c(float** c,int x_len,int y_len){
  float** res = new float*[y_len];
  for(int i =  0; i<y_len; i++ ){
    res[i] = new float[x_len];}
  float** f = get_matrix();
  for(int i =0; i <y_len; i++ ){
    for(int j = 0; j<x_len; j++){
      res[i][j]=f[i][j]*c[i][j];}
    }
    return res;
}

Point_<int> Picture::center_of_pressure(){

  float** min_distance_to_others;
  int nb_pts_unchanged=0;
  float **distance_to_others = new float*[x_length];
  for ( int i = 0 ; i < x_length ; i ++ ) {
    distance_to_others[i] = new float[y_length];
  }
  for(int i = 0 ; i < x_length ; i++){
    for(int j = 0; j< y_length ; j++){
      distance_to_others[i][j]=0;
      nb_pts_unchanged++;
      if(get_intensity(j,i)>=1){
        set_intensity(j,i,1);
      }
    }
  }
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      for (int ii=0;ii<x_length;ii++){
        for (int jj=0;jj<y_length;jj++){
          if(float p=get_intensity(jj,ii)!=1 && (ii!=i && jj!=j) ){
            distance_to_others[i][j]=distance_to_others[i][j]+sqrt(abs(ii-i)^2+abs(jj-j)^2)*(1-p);

          }
        }
      }
      distance_to_others[i][j]=distance_to_others[i][j]/(nb_pts_unchanged-1) ;
    }
  }
  int x_coordinate;
  int y_coordinate;
  float min_distance=distance_to_others[0][0];
  for (int i=0;i<x_length;i++){
    for (int j=0;j<y_length;j++){
      if(min_distance>distance_to_others[i][j]){
        x_coordinate=i;
        y_coordinate=j;
        min_distance=distance_to_others[i][j];
      }
    }
  }
Point p(x_coordinate,y_coordinate);
return(p);
}
