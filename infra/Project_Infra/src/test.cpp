#include "Picture.h"
#include "Usefull_functions.h"
#include <math.h>
#include<opencv2/opencv.hpp>
#include<iostream>
using namespace std;
using namespace cv;

float** ellipse_neighborhood(float** f, double a, double b, int x, int y, int x_len,int y_len){
  float** res = new float*[y_len];
  for(int i =  0; i<y_len; i++ ){
    res[i] = new float[x_len];}
  for(int i =0; i <y_len; i++ ){
    for(int j = 0; j<x_len; j++){
      if ((double((i-x)^2)/(pow(a,2))+(double((j-y)^2))/(pow(b,2)))<=1) {
        res[i][j]=f[i][j];
      }
      }
    }
    return res;

}


float** exp_c(int x_len, int y_len){
  float** res = new float*[y_len];
  for(int i =  0; i<y_len; i++ ){
    res[i] = new float[x_len];}
    for(int i =0; i <y_len; i++ ){
      for(int j = 0; j<x_len; j++){
        res[i][j]=exp(-sqrt(i^2+j^2));
      }}
  return res;
}

int main(int argc, char const *argv[]) {
  Mat image2 = imread( "../../data/strong_finger.png" , IMREAD_GRAYSCALE );

  Picture img2(image2);
  Point h=  img2.center_of_pressure();
  cout << h;
  img2.print_picture();
  for (int i=h.x-10;i<h.x+10;i++){
    for (int j= h.y -10; j< h.y+10; j++){
      img2.set_intensity(j,i,0.5);
    }
  }

  img2.print_picture();
  //img2.set_intensity(h.x,h.y,1);
  //cout<< "le point:"<< h << "\n";
  // float** f = img2.get_matrix();
  // display_matrix(f,img2.get_y_len(),img2.get_x_len());
  // float** g = ellipse_neighborhood(f,1,1,25,25,img2.get_x_len(),img2.get_y_len());
  // display_matrix(g,img2.get_y_len(),img2.get_y_len());
  // display_matrix(g,img2.get_y_len(),img2.get_x_len());
  // float** c = exp_c(img2.get_y_len(), img2.get_x_len());
  // float** res = img2.apply_c(c,img2.get_y_len(),  img2.get_x_len());
  // display_matrix(res, img2.get_y_len(), img2.get_x_len());
  return 0;
}
