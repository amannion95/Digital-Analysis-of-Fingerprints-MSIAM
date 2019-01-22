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

Point rotation_ij(Point a ,float rotation){
  a.x=a.x*cos(rotation)-a.y*sin(rotation);
  a.y=a.x*sin(rotation)+a.y*cos(rotation);
  cout << " a = " << a;
  return a;
}

float fct_c_test(float x){
  return 300./(300.+x);
  //return 1/(log(exp(1)+0.0000001*x));
  //return pow(1.000000001,-x);
}

vector<Point> segment(Point p0, Point p1){
  int dx = p1.x-p0.x;
  int dy = p1.y-p0.y;
  int nx = abs(dx);
  int ny = abs(dy);
  int sign_x = dx > 0? 1 : -1;
  int sign_y = dy > 0? 1 : -1;
  Point p(p0.x,p0.y);
  vector<Point> line;
  line.push_back(p);
  //cout << "nx = " << nx << " ny = " << ny;
  for(int ix =0, iy=0; ix<nx || iy < ny;){
    if((0.5+ix)/nx <= (0.5+iy)/ny){
      // the step is horizontal
      p.x +=sign_x;
      ix++;
    }
    else{
      // the step is vertical
      p.y+=sign_y;
      iy++;
    }
    line.push_back(Point(int(p.x),int(p.y)));
  }
  return line;
}

bool compare_polar_cord(Point2f a, Point2f b){
  return(a.y<b.y);
}

bool compare_y_cord(Point a, Point b){
  if(a.y != b.y){
    return(a.y<b.y);
  }
  else{
    return(a.x<b.x);
  }
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

Point2d pt_polar_rotation(Point2d orig, Point2d centre, double angle){
  if(orig==centre){
    return orig;
  }
  else{
    double rad = angle * M_PI / 180;
    double x = orig.x - centre.x;
    double y = centre.y - orig.y;

    double r = norm(Point(x,y));
    double theta = atan2(y,x);

    theta += rad;

    double rx = r * cos(theta) + (double)centre.x;
    double ry = (double)centre.y - r * sin(theta);

    return Point2d(rx,ry);
  }
}
