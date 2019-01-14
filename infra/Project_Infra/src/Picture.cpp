#include "Picture.h"
#include "Usefull_functions.h"
#include <vector>
#include <opencv2/opencv.hpp>
#include <math.h>
#include <cstdlib> // absolute value
#include <vector>

using namespace cv;
using namespace std;

Picture::Picture(const std::string& filename){
  picture=imread(filename,  IMREAD_GRAYSCALE);
  x_length=(picture.size()).width;
  y_length=(picture.size()).height;

}

Picture::Picture(unsigned int x_length,unsigned int y_length){
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
    std::cerr<<"Wrong intensity value, it must belong to [0,1]"<<std::endl;
  }
  picture.at<uchar>(i,j)=iftoi(intensity);
}

void Picture::print_picture()const{
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
  for(int i = 0 ; i < row ; i++ ){
    for(int j = 0 ; j < col ; j++ ){
      matrix[i][j] = iitof(picture.data[i*row+j]);
    }
  }
  return matrix;
}

Point_<int> Picture::center_of_pressure(){
  float threshold=0.1;
  Picture pressure_pic = this->clone();
  vector<Point> point_threshold;
  for(int i = 0 ; i < x_length ; i++){
    for(int j = 0; j< y_length ; j++){
      if(pressure_pic.get_intensity(j,i)>=threshold){
        pressure_pic.set_intensity(j,i,1);
      }
      else{
        //cout << get_intensity(j,i) << "    " ;
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
    else if(min(min_distance,distance[i])==distance[i]){
      indice=i;
      min_distance=distance[i];
    }
  }
return(point_threshold[indice]);
}

//-----------------------------------------------------------
//add by tristan 10th jan

//win_size must be odd ! Apply gaussian filter to the picture
Picture Picture::apply_gaussian_blur(int win_size=5)const{
  Mat blured_picture;
  GaussianBlur(picture,blured_picture,Size(win_size,win_size),0,0);
  Picture blured_Pic(blured_picture);
  return blured_Pic;
}

Point Picture::get_index_maximum_intensity()const{
  Point coord_min(0,0);
  Point coord_max(0,0);
  double x,y;
  minMaxLoc(picture, &x, &y, &coord_min, &coord_max);
  return coord_max;
}

Point Picture::get_index_minimum_intensity()const{
  Point coord_min(0,0);
  Point coord_max(0,0);
  double x,y;
  minMaxLoc(picture, &x, &y, &coord_min, &coord_max);
  return coord_min;
}

void Picture::print_pression_center(int size_win_gauss=5)const{
  int x,y;
  Picture img=apply_gaussian_blur(size_win_gauss);
  Picture print(picture);
  img.print_picture();
  x=img.get_index_minimum_intensity().x;
  y=img.get_index_minimum_intensity().y;
  for (int i=x-10;i<x+10;i++){
    for(int j=y-10;j<y+10;j++){
      img.set_intensity(j,i,0.7);
      print.set_intensity(j,i,0.7);
    }
  }
  img.print_picture();
  print.print_picture();
}

Point Picture::pressure_center_gauss(){
  Picture img;
  img=apply_gaussian_blur(31);
  return img.get_index_minimum_intensity();
}

//----------------------------------------------------------------------------

vector<Point> Picture::ellipse_nbh(Point p, unsigned int a, unsigned int b){
  int c1 = p.x;
  int c2 = p.y;
  vector<Point> nbh;
  for(unsigned int i=c1-a; i<=c1+a; i++){
    for(unsigned int j=c2-b; j<=c2+b; j++){
      if((double)((i-c1)*(i-c1)/pow(a,2) + (j-c2)*(j-c2)/pow(b,2)) <= 1){
        nbh.push_back(Point(i,j));
      }
    }
  }
  return nbh;
}

void Picture::show_nbh(vector<Point> nbh)const{
  Picture pic_w_nbh;
  pic_w_nbh = this->clone();
  for(Point &p : nbh){
    pic_w_nbh.set_intensity(p.y, p.x, 1);
  }
  pic_w_nbh.print_picture();
}

Picture Picture::log_transform_isotropic(Point p, unsigned int a, unsigned int b, double coef){
  Picture pressure_pic = this->clone();
  vector<Point> nbh = ellipse_nbh(p, a, b);
  for(Point &pixel : nbh){
    float c = log_coeff_isotropic(pixel, p, coef);
    float m = pressure_pic.get_intensity(pixel.y, pixel.x);
    pressure_pic.set_intensity(pixel.y, pixel.x, c*m);
  }
  return pressure_pic;
}

Picture Picture::pow_transform_isotropic(Point p, unsigned int a, unsigned int b, double coef){
  Picture pressure_pic = this->clone();
  vector<Point> nbh = ellipse_nbh(p, a, b);
  for(Point &pixel : nbh){
    float r = norm(p-pixel)/(max(a,b));
    float m = pressure_pic.get_intensity(pixel.y, pixel.x);
    pressure_pic.set_intensity(pixel.y, pixel.x, sqrt(r)*m);
  }
  return pressure_pic;
}

//----------------------------------------------------


Picture Picture::extract_ellipse_pic(Point center, unsigned int a,unsigned int b){
  int c1=center.x;
  int c2=center.y;

  Picture extraction=clone();
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      if((double)((i-c1)*(i-c1)/pow(a,2) + (j-c2)*(j-c2)/pow(b,2)) > 1){
        extraction.set_intensity(j,i,1);
      }
    }
  }
  return extraction;
}


Picture Picture::apply_threshold(float set_lim){
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


vector<Point> Picture::get_0intensity_index (){
  vector<Point> contain;
  for(int i=0;i<x_length;i++){
    for (int j=0;j<y_length;j++){
      if (get_intensity(j,i)==0){
        Point p(j,i);
        contain.push_back(p);
      }
    }
  }
  return contain;
}

Point Picture::get_median_center(vector<Point> intensity_index){
  vector<int> x_index;
  vector<int> y_index;
  Point coord(0,0);

  for (vector<Point>::iterator i = intensity_index.begin() ; i != intensity_index.end(); i++){
    x_index.push_back(i->x);
    y_index.push_back(i->y);
  }


  sort(x_index.begin(),x_index.end());
  sort(y_index.begin(),y_index.end());

  if (intensity_index.size()==0){
    cerr<<"Empty vector do not have any median value."<<endl;
    return coord;
  }
  coord.x=x_index[intensity_index.size() / 2];
  coord.y=y_index[intensity_index.size() / 2];
  return coord;
}

//BAD RESULTS, NOT USABLE
void Picture::print_median_center(int thresh=0.01){
  int x,y;
  Picture img=apply_threshold(thresh);
  Picture print(picture);
  Point* p=new Point(img.get_median_center(img.get_0intensity_index()));
  x=p->x;
  y=p->y;

  for (int i=x-10;i<x+10;i++){
    for(int j=y-10;j<y+10;j++){
      img.set_intensity(j,i,0.7);
      print.set_intensity(j,i,0.7);
    }
  }
  delete p;
  img.print_picture();
  print.print_picture();

}


Picture Picture::try_apply_isotrop(cv::Point center,unsigned int a,unsigned int b){

  Picture extrac_smoothed;
  Point p(0,0);
  extrac_smoothed=extract_ellipse_pic(center,a,b);
  float distance;
  for(int i=0;i<extrac_smoothed.x_length;i++){
    for (int j=0;j<extrac_smoothed.y_length;j++){
      if (extrac_smoothed.get_intensity(j,i)!=1)
        p.y=j;
        p.x=i;
        distance=norm(center-p);
        extrac_smoothed.set_intensity(j,i,fct_c_test(distance)*extrac_smoothed.get_intensity(j,i));
    }
  }
  return extrac_smoothed.extract_ellipse_pic(center,a,b);
}
