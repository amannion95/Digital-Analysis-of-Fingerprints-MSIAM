#include "Useful_functions.h"
#include <iostream>
#include <math.h>
#include <opencv2/opencv.hpp>


using namespace std;
using namespace cv;
using namespace Eigen;

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
  if(a.y != b.y){
    return(a.y<b.y);
  }
  else{
    return(a.x<b.x);
  }
}


bool compare_y_cord(Point a, Point b){
  if(a.y != b.y){
    return(a.y<b.y);
  }
  else{
    return(a.x<b.x);
  }
}
void convolveDFT(InputArray A, InputArray B, OutputArray C)
{
    // reallocate the output array if needed
    C.create(abs(A.rows() + B.rows()) - 1, abs(A.cols() + B.cols()) - 1, A.type());
    Size dftSize;
    // calculate the size of DFT transform
    dftSize.width = getOptimalDFTSize(A.cols() + B.cols() - 1);
    dftSize.height = getOptimalDFTSize(A.rows() + B.rows() - 1);

    // allocate temporary buffers and initialize them with 0's
    Mat tempA(dftSize, A.type(), Scalar::all(0));
    Mat tempB(dftSize, B.type(), Scalar::all(0));
    // copy A and B to the top-left corners of tempA and tempB, respectively
    Mat roiA(tempA, Rect(0,0,A.cols(),A.rows()));
    A.copyTo(roiA);
    Mat roiB(tempB, Rect(0,0,B.cols(),B.rows()));
    B.copyTo(roiB);
    // now transform the padded A & B in-place;
    // use "nonzeroRows" hint for faster processing
    dft(tempA, tempA, 0, A.rows());
    dft(tempB, tempB, 0, B.rows());
    // multiply the spectrums;
    // the function handles packed spectrum representations well
    mulSpectrums(tempA, tempB, tempA,0,false);
    // transform the product back from the frequency domain.
    // Even though all the result rows will be non-zero,
    // you need only the first C.rows of them, and thus you
    // pass nonzeroRows == C.rows
    dft(tempA, tempA, DFT_INVERSE + DFT_SCALE, C.rows());
    // now copy the result back to C.
    tempA(Rect(0, 0, C.cols(), C.rows())).copyTo(C);
    // all the temporary buffers will be deallocated automatically
}

Matrix<float,Dynamic,Dynamic> Evolutive_kernel_no_identity(float t, int size){

  Matrix<float,Dynamic,Dynamic> Identity(size,size);
  Matrix<float,Dynamic,Dynamic> mask(size,size);
  Matrix<float,Dynamic,Dynamic> evolutive_mask(size,size);
  Point center(mask.rows()/2,mask.cols()/2);
  float msqe=pow((1-t),1);

  for(int i = 0 ; i<size;i++){
    for(int j=0 ; j< size;j++){
      Identity(j,i)=0.;
      if(i==size/2 && j==size/2){
        Identity(j,i)=1.;
      }
      mask(j,i)=(1/(2*3.141592))*exp(-(pow(j-center.y,2)+pow(i-center.x,2))/2*msqe);
      evolutive_mask(j,i)=mask(j,i);
    }
  }


  return(evolutive_mask);
}

Matrix<float,Dynamic,Dynamic> Evolutive_kernel(float t, int size){

  Matrix<float,Dynamic,Dynamic> Identity(size,size);
  Matrix<float,Dynamic,Dynamic> mask(size,size);
  Matrix<float,Dynamic,Dynamic> evolutive_mask(size,size);
  Point center(mask.rows()/2,mask.cols()/2);
  float msqe=pow((1-t),1);

  for(int i = 0 ; i<size;i++){
    for(int j=0 ; j< size;j++){
      Identity(j,i)=0.;
      if(i==size/2 && j==size/2){
        Identity(j,i)=1.;
      }
      mask(j,i)=(1/(2*3.141592*pow(msqe,1.5)))*exp(-(pow(j-center.y,2)+pow(i-center.x,2))/2*msqe);
      evolutive_mask(j,i)=(1-pow(t,1))*Identity(j,i)+pow(t,1)*mask(j,i);
    }
  }


  return(evolutive_mask);
}

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
