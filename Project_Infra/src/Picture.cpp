#include "Picture.h"


using namespace cv;
using namespace std;
using namespace Eigen;

Picture::Picture(const std::string& filename){
  picture=imread(filename,  IMREAD_GRAYSCALE);
  x_length=(picture.size()).width;
  y_length=(picture.size()).height;
}

Picture::Picture(unsigned int x_length,unsigned int y_length){
  Mat image(y_length,x_length,CV_8UC1,Scalar(255));
  picture=image.clone();
  this->x_length=x_length;
  this->y_length=y_length;
}

Picture::Picture(const cv::Mat& pic){
  picture=pic.clone();
  x_length=(picture.size()).width;
  y_length=(picture.size()).height;
}

Picture::Picture(){
  x_length=0;
  y_length=0;
  Mat image(0,0,CV_8UC1,Scalar(255));
  picture=image;
}

float Picture::get_intensity(unsigned int i, unsigned int j)const{
  return iitof((int)picture.at<uchar>(i,j));
}

void Picture::set_intensity(unsigned int i, unsigned int j,float intensity){
  if ((intensity<0)||(intensity>1)){
    std::cerr<<"Wrong intensity value, it must belong to [0,1]"<<std::endl;
  }
  picture.at<uchar>(i,j)=iftoi(intensity);
}

void Picture::print_picture(string name)const{
  namedWindow(name, WINDOW_NORMAL );
  imshow(name, picture);
  waitKey(0);
}

void Picture::SAVE_PIC(string name)const{
  imwrite(name, picture);
}

void Picture::operator=(Picture Pic){
  picture=(Pic.picture).clone();
  x_length=Pic.x_length;
  y_length=Pic.y_length;
}

Picture Picture::operator-(Picture to_substract)const{
  float ab;
  Picture diff_pic(x_length, y_length);
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      ab=abs(get_intensity(j,i)-to_substract.get_intensity(j,i));
      diff_pic.set_intensity(j,i,ab);
    }
  }
  return diff_pic;
}

float Picture::maximum_intensity()const{
  float max_intensity=0.;
  for (int i = 0 ; i < x_length ; i++ ){
    for (int j = 0 ; j < y_length ; j++ ){
      if (max_intensity<get_intensity(j,i)){
          max_intensity=get_intensity(j,i);
      }
    }
  }
  return(max_intensity);
}

float Picture::minimum_intensity()const{
  float min_intensity=1.;
  for (int i = 0 ; i < x_length ; i++ ){
    for (int j = 0 ; j < y_length ; j++ ){
      if (min_intensity>get_intensity(j,i)){
          min_intensity=get_intensity(j,i);
      }
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

unsigned int Picture::get_x_len()const{
  return x_length;
}
unsigned int Picture::get_y_len()const{
  return y_length;
}

void Picture::rescale_color(){
  float min_intensity = minimum_intensity();
  float max_intensity = maximum_intensity();
  float size = max_intensity - min_intensity;
  assert ( size != 0 ) ;
  for (int i = 0 ; i < x_length ; i++ ){
    for (int j = 0 ; j < y_length ; j++ ){
      set_intensity(j,i,(get_intensity(j,i)-min_intensity)/size);
    }
  }
}
Matrix<float,Dynamic,Dynamic> Picture::get_matrix()const{
  int col=x_length;  //= 256
  int row=y_length;  //= 288
  Matrix<float,Dynamic,Dynamic> mat(row,col);
  for(int i = 0 ; i < col; i++ ){
    for(int j = 0 ; j < row ; j++ ){
      mat(j,i) = iitof(picture.data[i*row+j]);
    }
  }
  return mat;
}

Point Picture::center_of_pressure()const{
  float threshold=0.1;
  Picture pressure_pic = clone();
  vector<Point> point_threshold;
  for(int i = 0 ; i < x_length ; i++){
    for(int j = 0; j< y_length ; j++){
      if(pressure_pic.get_intensity(j,i)>=threshold){
        pressure_pic.set_intensity(j,i,1);
      }
      else{
        point_threshold.push_back(Point(i,j));
      }
    }
  }
  vector<float> distance;
  float min_distance=0;
  int indice;
  for(int i=0;i<point_threshold.size();i++){
    distance.push_back(0);
    for(int j=0;j<point_threshold.size();j++){
      distance[i]=distance[i]+float(norm(point_threshold[i]-point_threshold[j]));
    }
    if(i==0){
      min_distance=distance[i];
    }
    else if(distance[i]<min_distance){
      indice=i;
      min_distance=distance[i];
    }
  }
  return(point_threshold[indice]);
}

//win_size must be odd ! Apply gaussian filter to the picture
Picture Picture::apply_gaussian_blur(int win_size=5)const{
  Mat blured_picture;
  GaussianBlur(picture,blured_picture,Size(win_size,win_size),0,0);
  Picture blured_Pic(blured_picture);
  return blured_Pic;
}

Point Picture::get_index_minimum_intensity()const{
  Point coord_min(0,0);
  Point coord_max(0,0);
  double x,y;
  minMaxLoc(picture, &x, &y, &coord_min, &coord_max);
  return coord_min;
}

void Picture::print_pression_center_gauss_threshold(){
  int x,y;

  Picture img=apply_threshold(0.1).apply_gaussian_blur(51);
  Picture print(picture);
  x=img.get_index_minimum_intensity().x;
  y=img.get_index_minimum_intensity().y;
  for (int i=x-10;i<x+10;i++){
    for(int j=y-10;j<y+10;j++){
      print.set_intensity(j,i,0.7);
    }
  }
  print.print_picture();
}

Point Picture::pressure_center_gauss_threshold()const{
  Picture img=apply_threshold(0.1);
  img=img.apply_gaussian_blur(51);
  return img.get_index_minimum_intensity();
}

Picture Picture::apply_threshold(float set_lim)const{
  Picture th_ed = clone();
  for(int i =0; i<x_length; i++){
    for(int j =0; j<y_length; j++){
      if(get_intensity(j,i)>set_lim){
        th_ed.set_intensity(j,i,1);
      }
      else{
        th_ed.set_intensity(j,i,0);
      }
    }
  }
  return th_ed;
}

bool Picture::isinframe(Point p)const{
  return p.x>=0 && p.x<x_length &&p.y>=0 && p.y<y_length;
}
