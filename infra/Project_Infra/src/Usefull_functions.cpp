#include "Usefull_functions.h"

#include <math.h>
#include <opencv2/opencv.hpp>
#include <iostream>


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
    return(-1);}
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
  return c * log(r+1);
}

float pow_coeff_isotropic(Point centre, Point p, int n, double c){
  float r = norm(p-centre);
  return c / pow(r+1, n);
}


float symmetry_wrt_05(float x){
  if (x<0.5){
    return x+2*abs(x-0.5);
  }
  else{
    return x-2*abs(x-0.5);
  }
}



float fct_c_test(float x){
  return 500./(500.+x);
}
