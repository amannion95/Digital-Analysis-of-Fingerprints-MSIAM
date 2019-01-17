#include <iostream>
#include <cstdlib>
#include <math.h>
#include <opencv2/opencv.hpp>
#include <vector>
#include <algorithm>

#include "Picture.h"
#include "Usefull_functions.h"

using namespace std;
using namespace cv;

/*this function creates rotated pixels by converting to polar co-ordinates and
adding to the angle, then assigns pixel intensities to the target frame by simply
casting the rotated pixel co-ordinates to integers*/
Picture Picture::cast_rotation_polar(Point centre, double angle){
  Picture rotated_pic(x_length, y_length);
  double rad = -angle * M_PI/180;
  for(int i=0; i<x_length; i++){
    for(int j=0; j<y_length; j++){
      int x = i - centre.x;
      int y = centre.y - j;

      //polars
      double r = norm(Point(x,y));
      double theta = atan2((double)y, (double)x);

      theta += rad;

      //convert back to original system
      int rx = (int)(round(r * cos(theta))) + centre.x;
      int ry = centre.y - (int)(round(r * sin(theta)));
      //cout << Point(i,j) << "   " << Point(rx, ry) << endl;

      if(rx >= 0 && rx < x_length && ry >= 0 && ry < y_length){
        rotated_pic.set_intensity(j, i, get_intensity(ry, rx));
      }
      else{
        rotated_pic.set_intensity(j, i, 1.0);
      }
    }
  }
  return rotated_pic;
}

/*Rotation using nearest-neighbour interpolation to assign intensities; finds the
minimum
*/
Picture Picture::nn_rotation_polar(Point centre, double angle){
  Picture rotated_pic(x_length, y_length);
  double rad = -angle * M_PI/180;
  for(int i=0; i<x_length; i++){
    for(int j=0; j<y_length; j++){
      int x_ = i - centre.x;
      int y_ = centre.y - j;

      //polars
      double r = norm(Point(x_,y_));
      double theta = atan2((double)y_, (double)x_);

      theta += rad;

      //convert back to original system
      double x = r * cos(theta) + (double)(centre.x);
      double y = (double)(centre.y) - r * sin(theta);
      //cout << Point(i,j) << "   " << Point(x, y) << endl;

      if(x >= 0 && x < x_length && y >= 0 && y < y_length){
        Point2d target(x,y);
        double xf = floor(x);
        double xc = ceil(x);
        double yf = floor(y);
        double yc = ceil(y);

        //calculate Euclidean distances, check for minimum and assign intensity
        //at the corresponding point in the initial image
        double d1 = norm(target-Point2d(xf,yf));
        double d2 = norm(target-Point2d(xf,yc));
        double d3 = norm(target-Point2d(xc,yf));
        double d4 = norm(target-Point2d(xc,yc));

        double distances[] = {d1, d2, d3, d4};
        double min_dist = *min_element(distances, distances+4);

        if(min_dist==d1){
          rotated_pic.set_intensity(j, i, get_intensity(yf,xf));
        }
        else if(min_dist==d2){
          rotated_pic.set_intensity(j, i, get_intensity(yc,xf));
        }
        else if(min_dist==d3){
          rotated_pic.set_intensity(j, i, get_intensity(yf,xc));
        }
        else{
          rotated_pic.set_intensity(j, i, get_intensity(yc,xc));
        }
      }
      else{
        rotated_pic.set_intensity(j, i, 1.0);
      }
    }
  }
  return rotated_pic;
}

Picture Picture::bilinear_rotation_polar(Point centre, double angle){
  Picture rotated_pic(x_length, y_length);
  double rad = -angle * M_PI/180;
  for(int i=0; i<x_length; i++){
    for(int j=0; j<y_length; j++){
      int x_ = i - centre.x;
      int y_ = centre.y - j;

      //polars
      double r = norm(Point(x_,y_));
      double theta = atan2((double)y_, (double)x_);

      theta += rad;

      //convert back to original system
      double x = r * cos(theta) + (double)(centre.x);
      double y = (double)(centre.y) - r * sin(theta);
      //cout << Point(i,j) << "   " << Point(x, y) << endl;

      if(x >= 0 && x < x_length && y >= 0 && y < y_length){
        int xf = (int)(floor(x));
        int xc = (int)(ceil(x));
        int yf = (int)(floor(y));
        int yc = (int)(ceil(y));

        if(xf==xc || yf==yc){
          rotated_pic.set_intensity(j, i, get_intensity((int)y, (int)x));
        }

        else{
          //step 1: interpolation in x-direction
          float lower_left = get_intensity(yf, xf);
          float lower_right = get_intensity(yf, xc);
          float i1 = (xc - x)*lower_left + (x - xf)*lower_right;

          float upper_left = get_intensity(yc, xf);
          float upper_right = get_intensity(yf, xf);
          float i2 = (xc - x)*upper_left + (x - xf)*upper_right;

          //step 2: interpolation in y-direction
          float I = (yc - y)*i1 + (y - yf)*i2;

          rotated_pic.set_intensity(j, i, I);
          }
        }
      else{
        rotated_pic.set_intensity(j, i, 1.0);
      }
    }
  }
  return rotated_pic;
}

vector<Point2d> rotated_pixels(Point centre, double angle){
  vector<Point2d> co_ords;
  for(unsigned int i = 0; i < get_x_len(); i++){
    for(unsigned int j = 0; j < get_y_len(); j++){
      Point2d r = rotation(Point(i,j), centre, angle);
      co_ords.push_back(r);
    }
  }
  cout << "pixels rotated" << endl;
  return co_ords;
}
