#include "Usefull_functions.h"

#include <math.h>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <cassert>

using namespace std;
using namespace cv;

float iitof(int color){
  if(color<0 || 255<color){
    cout<<"invalid color"<<"\n";
    return(-1);}
    else return(float(color)/255);
}

int iftoi(float color){
  if(color<0 || 1<color){
    cout<<"invalid color"<<"\n";
    return(-1);
  }
    else return(int(color*255));
}

void display_matrix(float ** matrix,int row, int col){
  for (int i=0; i<col;i++){
    for(int j=0; j<row;j++){
      cout << matrix[i][j] << "\t";
    }
    cout << "\n";
  }
}

float log_coeff_isotropic(Point centre, Point p, double c){
  float r = norm(p-centre);
  return c * log(1+r);
}

/*float pow_coeff_isotropic(Point centre, Point p, int n, float c){
  float r = norm(p-centre);
  return c * (pow(r,n) + 1);
}

float intensity_sym(float i){
  assert(i <= 1);
  if((i-0.5) < 0){
    return i + 2 * abs(i-0.5);
  }
  else{
    return i - 2 * abs(i-0.5);
  }
}*/
