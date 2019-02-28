#include <iostream>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <math.h>
#include <vector>
#include <algorithm>
#include <cassert>

#include "Picture.h"

using namespace std;
using namespace cv;
using namespace Eigen;

Picture Picture::cast_rotation_cart(Point centre, double angle)const{
  Picture rotated_pic(x_length, y_length);
  double rad = M_PI/(180)*angle;
  for(int i=0; i<x_length; i++){
    for(int j=0; j<y_length; j++){
      int x = int(cos(rad)*(i-centre.x) -sin(rad)*(j-centre.y)) + centre.x;
      int y = int(sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y)) + centre.y;
      if(isinframe(Point(x,y))){
        rotated_pic.set_intensity(j, i, get_intensity(y,x));
      }
    }
  }
  return rotated_pic;
}

Picture Picture::cast_rotation_polar(Point centre, double angle)const{
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

      if(isinframe(Point(rx,ry))){
        rotated_pic.set_intensity(j, i, get_intensity(ry, rx));
      }
      else{
        rotated_pic.set_intensity(j, i, 1.0);
      }
    }
  }
  return rotated_pic;
}

Picture Picture::nn_rotation_polar(Point& centre, double angle)const{
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

Picture Picture::bilinear_rotation_polar(Point& centre, double angle)const{
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

      if(isinframe(Point(x,y))){
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
          float upper_right = get_intensity(yc, xc);
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

float Picture::bilinear_interpolation(Point2d p)const{
  if(isinframe(p)){
    int xf = (int)floor(p.x);
    int xc = (int)ceil(p.x);
    int yf = (int)floor(p.y);
    int yc = (int)ceil(p.y);

    if(xf==xc || yf == yc){
      return get_intensity((int)p.y, (int)p.x);
    }
    else{
      float i1 = (xc - p.x) * get_intensity(yf,xf)
                  + (p.x - xf) * get_intensity(yf,xc);
      float i2 = (xc - p.x) * get_intensity(yc,xf)
                  + (p.x - xf) * get_intensity(yc,xc);
      float I = (yc - p.y) * i1 + (p.y - yf) * i2;
      return I;
    }
  }
  else{
    return 0;
  }
}

Point Picture::local_cop(list<Point>& region)const{
  float threshold = 0.1;
  Picture pp = clone();
  vector<Point> threshpoints;
  for(Point& p : region){
    if(pp.get_intensity(p.y, p.x) >= threshold){
      pp.set_intensity(p.y, p.x, 1);
    }
    else{
      threshpoints.push_back(p);
    }
  }
  vector<float> dist_sums;
  for(Point& p : threshpoints){
    float dsum = 0.0;
    for(Point& q : threshpoints){
      dsum += norm(p-q);
    }
    dist_sums.push_back(dsum);
  }
  vector<float>::iterator min_dist = min_element(dist_sums.begin(), dist_sums.end());
  int index = distance(dist_sums.begin(), min_dist);
  return threshpoints[index];
}

Matrix<list<Point>, Dynamic, Dynamic> Picture::morph_subsets(int radius, int subset_sidelen, Point centre)const{
  Matrix<list<Point>, Dynamic, Dynamic> subsets(subset_sidelen, subset_sidelen);
  int zero_x = centre.x - radius;
  int zero_y = centre.y - radius;
  cout << "zero: " << Point(zero_x, zero_y) << endl;
  int subset_size = 2*radius/subset_sidelen;
  cout << "subset size: " << subset_size <<endl;
  for(size_t i=0; i<subset_sidelen; i++){
    for(size_t j=0; j<subset_sidelen; j++){
      //cout << "checkpoint 1" << endl;
      list<Point> subset;
      //cout << "x: " << zero_x+i*subset_size << endl;
      for(unsigned int k=zero_y+i*subset_size;
          k<=zero_y+i*subset_size+subset_size; k++){
            for(unsigned int l=zero_x+j*subset_size;
                l<=zero_x+j*subset_size+subset_size; l++){
                  subset.emplace_back(Point(j+l, i+k));
            }
      }
      subsets(i,j) = subset;
    }
  }
  cout << "subsets formed" << endl;
  return subsets;
}

double Picture::local_cop_distance(list<Point> subset, Point centre)const{
  Point cop = local_cop(subset);
  double dist = norm(cop - centre);
  return dist;
}

Picture Picture::swirl(Point centre, double twist, int radius)const{
  Picture swirl_pic(x_length, y_length);
  for(unsigned int i=0; i<x_length; i++){
    for(unsigned int j=0; j<y_length; j++){
      double dist = norm(centre - Point(i,j));
      if(dist < radius){
        double amt = 1 - dist / (double)radius;
        double angle = 2 * twist * M_PI * amt;
        Point2d p = pt_polar_rotation(Point2d((double)i, (double)j), centre, -angle);
        float intensity = this->bilinear_interpolation(p);
        swirl_pic.set_intensity(j, i, intensity);
      }
      else{
        swirl_pic.set_intensity(j, i, get_intensity(j, i));
      }
    }
  }
  return swirl_pic;
}

Picture Picture::local_erosion(Point centre, int a, int b)const{
  Picture eroded_pic = clone();
  vector<Point> domain = ellipse_nbh(centre, a, b);
  //Picture temp_pic = clone();
  for(Point& p : domain){
    if(isinframe(p)){
      Vector3f x_coef, y_coef;
      float dx = (float)abs(p.x-centre.x);
      float dy = (float)abs(p.y-centre.y);
      for(int i=-1; i<2; i++){
        float d1 = (dx+i)/(float)a;
        float d2 = (dy+i)/(float)b;
        x_coef(i+1) = 1/(sqrt(2*M_PI))*(exp(-dx*dx*0.5)-1);
        y_coef(i+1) = 1/(sqrt(2*M_PI))*(exp(-dy*dy*0.5)-1);
      }
      Matrix3f structelem = x_coef * y_coef.transpose();
      //cout << "structuring element: \n" << structelem << endl;
      float energy = structelem.sum();
      if(energy != 0){
        structelem /= energy;
      }
      Matrix3f eroded_intensities;
      float fa = (float)a;
      float fb = (float)b;
      float d = norm(Point2f(dx/fa, dy/fb));
      for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
          eroded_intensities(i,j) = get_intensity(p.y-1+j, p.x-1+i) - structelem(i,j);
        }
      }
      float chosen_intensity = d*d*get_intensity(p.y, p.x) + (1.0-d*d)*eroded_intensities.minCoeff();
      if(chosen_intensity < 0.0){
        chosen_intensity = 0.0;
      }
      eroded_pic.set_intensity(p.y, p.x, chosen_intensity);
    }
  }
  return eroded_pic;
}

Picture Picture::local_dilation(Point centre, int a, int b)const{
  Picture dilated_pic = clone();
  vector<Point> domain = ellipse_nbh(centre, a, b);
  for(Point& p: domain){
    if(isinframe(p)){
      Vector3f x_coef, y_coef;
      float dx = (float)abs(p.x-centre.x);
      float dy = (float)abs(p.y-centre.y);
      for(int i=-1; i<2; i++){
        float dx = (dx+i)/(float)a;
        float dy = (dy+i)/(float)b;
        x_coef(i+1) = 1/(sqrt(2*M_PI))*(exp(-dx*dx*0.5)-1);
        y_coef(i+1) = 1/(sqrt(2*M_PI))*(exp(-dy*dy*0.5)-1);
      }
      Matrix3f structelem = x_coef * y_coef.transpose();
      float energy = structelem.sum();
      if(energy != 0){
        structelem /= energy;
      }
      Matrix3f dilated_intensities;
      float d = norm(Point2f(dx/(float)a, dy/(float)b));
      for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
          dilated_intensities(i,j) = get_intensity(p.y-1+j, p.x-1+i) + structelem(i,j);
        }
      }
      float chosen_intensity = d*d*get_intensity(p.y, p.x) + (1.0-d*d)*dilated_intensities.maxCoeff();
      if(chosen_intensity > 1.0){
        chosen_intensity = 1.0;
      }
      dilated_pic.set_intensity(p.y, p.x, chosen_intensity);
    }
  }
  return dilated_pic;
}

Picture Picture::swirl_zonalmorph(Point centre, double twist, int radius, int root_num_sub)const{
  Picture swirl_pic = swirl(centre, twist, radius);

  //construct matrix of subsets
  Matrix<list<Point>, Dynamic, Dynamic> subsets(root_num_sub, root_num_sub);
  int zero_x = centre.x - radius;
  int zero_y = centre.y - radius;
  int subset_size = 2*radius/root_num_sub;
  for(size_t i=0; i<root_num_sub; i++){
    for(size_t j=0; j<root_num_sub; j++){
      list<Point> subset;
      for(unsigned int k=zero_y+i*subset_size;
        k<=zero_y+(i+1)*subset_size; k++){
          for(unsigned int l=zero_x+j*subset_size;
            l<=zero_x+(j+1)*subset_size; l++){
              subset.emplace_back(Point(j+l, i+k));
          }
      }
      subsets(i,j) = subset;
    }
  }

  //construct "morphology indicator" matrix
  MatrixXi morph_indic(root_num_sub, root_num_sub);
  for(size_t j=0; j<subsets.cols(); j++){
    for(size_t i=0; i<subsets.rows(); i++){
      list<Point> subset = subsets(i,j);
      Point topl = subset.front();
      Point botr = subset.back();
      if(isinframe(topl) && isinframe(botr)){
        float d1 = this->local_cop_distance(subset, centre);
        float d2 = swirl_pic.local_cop_distance(subset, centre);
        if(d1 > d2){
          morph_indic(i,j) = 1;
        }
        else if(d1 < d2){
          morph_indic(i,j) = -1;
        }
        else{
          morph_indic(i,j) = 0;
        }
      }
      else{
        morph_indic(i,j) = 0;
      }
    }
  }
  cout << "morphology indicator: \n" << morph_indic << endl;

  //padded indicator matrix to avoid boundary awkwardness
  MatrixXi pad_indic(morph_indic.rows()+2, morph_indic.cols()+2);
  for(size_t l=0; l<pad_indic.cols(); l++){
    pad_indic(0,l) = 0;
    pad_indic(pad_indic.rows()-1,l) = 0;
  }
  for(size_t k=0; k<pad_indic.rows(); k++){
    pad_indic(k,0) = 0;
    pad_indic(k,pad_indic.cols()-1) = 0;
  }
  for(size_t l=1; l<pad_indic.cols()-1; l++){
    for(size_t k=1; k<pad_indic.rows()-1; k++){
      pad_indic(k,l) = morph_indic(k-1,l-1);
    }
  }

  //loop over each grid square...
  for(size_t j=0; j<subsets.cols(); j++){
    for(size_t i=0; i<subsets.rows(); i++){

      //find centre of square
      list<int> x_list, y_list;
      for(Point& p : subsets(i,j)){
        x_list.emplace_back(p.x);
        y_list.emplace_back(p.y);
      }
      int xmax = *max_element(x_list.begin(), x_list.end());
      int ymax = *max_element(y_list.begin(), y_list.end());
      int xmin = *min_element(x_list.begin(), x_list.end());
      int ymin = *min_element(y_list.begin(), y_list.end());
      int halflen = (int)round((xmax-xmin)/2);
      Point sub_cen(xmin+halflen, ymin+halflen);

      //construct boolean array to represent indicator neighbourhood; "trbl" stands
      //for "top, right, bottom, left"; the ordering of the 4-array. If the neighbour
      //above the point p has the same indicator value as p, the first entry in the array
      //will evaluate to false and so on.
      int neighbours_trbl[4] = {pad_indic(i,j+1), pad_indic(i+1,j+2),
                                pad_indic(i+2,j+1), pad_indic(i+1,j)};
      bool dec_trbl[4];
      for(unsigned int k=0; k<4; k++){
        if(neighbours_trbl[k] != morph_indic(i,j)){
          dec_trbl[k] = 1;
        }
        else{
          dec_trbl[k] = 0;
        }
      }

      //morphological filtering part
      Picture temp_pic = swirl_pic.clone();
      if(morph_indic(i,j) != 0){
        for(Point& p : subsets(i,j)){

          //cross-shape structuring element
          list<Point> cross = {p};
          for(int i=-1; i<2; i++){
            cross.emplace_back(Point(p.x+i, p.y));
            cross.emplace_back(Point(p.x, p.y+i));
          }
          list<float> structelem;
          for(Point& q : cross){
            structelem.emplace_back(temp_pic.get_intensity(q.y, q.x));
          }
          float chosen_intensity;
          if(morph_indic(i,j) == 1){
            chosen_intensity = *max_element(structelem.begin(), structelem.end());
          }
          else{
            chosen_intensity = *min_element(structelem.begin(), structelem.end());
          }
          assert(chosen_intensity>=0.0 && chosen_intensity<=1.0);

          float weighted_intensity = chosen_intensity;

          //polar co-ordinate angle with respect to the centre of the square
          float dx = (float)abs(p.x-sub_cen.x);
          float dy = (float)abs(p.y-sub_cen.y);
          double theta_p = atan2((double)dy, (double)dx);

          //diagonal angle marker
          double a = M_PI*0.25;

          //decreasing morphology effect weighting
          double to_edge;
          if(theta_p >= a && theta_p < 3.0*a){
            to_edge = norm(p-Point(p.x, ymin));
          }
          else if(theta_p >= 3.0*a && theta_p < 5.0*a){
            to_edge = norm(p-Point(xmin, p.y));
          }
          else if(theta_p >= 5.0*a && theta_p < 7.0*a){
            to_edge = norm(p-Point(p.x, ymax));
          }
          else{
            to_edge = norm(p-Point(xmax, p.y));
          }
          double d = norm(p-sub_cen)/(norm(p-sub_cen)+to_edge);
          float w = d*temp_pic.get_intensity(p.y, p.x) + (1.0-d)*chosen_intensity;

          //Sixteen possible neighbour configurations...
          if(dec_trbl[0] && dec_trbl[1] && dec_trbl[2] && dec_trbl[3]){
            weighted_intensity = w;
          }
          else if(!dec_trbl[0] && dec_trbl[1] && dec_trbl[2] && dec_trbl[3]){
            if(!(theta_p > a && theta_p < 3.0*a)){
              weighted_intensity = w;
            }
          }
          else if(dec_trbl[0] && !dec_trbl[1] && dec_trbl[2] && dec_trbl[3]){
            if(!(theta_p < a || theta_p > 7.0*a)){
              weighted_intensity = w;
            }
          }
          else if(dec_trbl[0] && dec_trbl[1] && !dec_trbl[2] && dec_trbl[3]){
            if(!(theta_p > 5.0*a && theta_p < 7.0*a)){
              weighted_intensity = w;
            }
          }
          else if(dec_trbl[0] && dec_trbl[1] && dec_trbl[2] && !dec_trbl[3]){
            if(!(theta_p > 3.0*a && theta_p < 5.0*a)){
              weighted_intensity = w;
            }
          }
          else if(!dec_trbl[0] && !dec_trbl[1] && dec_trbl[2] && dec_trbl[3]){
            if(theta_p > 3.0*a && theta_p < 5.0*a){
              weighted_intensity = w;
            }
          }
          else if(!dec_trbl[0] && dec_trbl[1] && !dec_trbl[2] && dec_trbl[3]){
            if(!((theta_p > a && theta_p < 3.0*a) || (theta_p > 5.0*a && theta_p < 7.0*a))){
              weighted_intensity = w;
            }
          }
          else if(!dec_trbl[0] && dec_trbl[1] && dec_trbl[2] && !dec_trbl[3]){
            if(!(theta_p > a && theta_p < 5.0*a)){
              weighted_intensity = w;
            }
          }
          else if(dec_trbl[0] && !dec_trbl[1] && !dec_trbl[2] && dec_trbl[3]){
            if(theta_p > a && theta_p < 5.0*a){
              weighted_intensity = w;
            }
          }
          else if(dec_trbl[0] && !dec_trbl[1] && dec_trbl[2] && !dec_trbl[3]){
            if((theta_p > a && theta_p < 3.0*a) || (theta_p > 5.0*a && theta_p < 7.0*a)){
              weighted_intensity = w;
            }
          }
          else if(dec_trbl[0] && dec_trbl[1] && !dec_trbl[2] && !dec_trbl[3]){
            if(!(theta_p > 3.0*a && theta_p < 7.0*a)){
              weighted_intensity = w;
            }
          }
          else if(!dec_trbl[0] && !dec_trbl[1] && !dec_trbl[2] && dec_trbl[3]){
            if(theta_p > 3.0*a && theta_p < 5.0*a){
              weighted_intensity = w;
            }
          }
          else if(!dec_trbl[0] && !dec_trbl[1] && dec_trbl[2] && !dec_trbl[3]){
            if(theta_p > 5.0*a && theta_p < 7.0*a){
              weighted_intensity = w;
            }
          }
          else if(!dec_trbl[0] && dec_trbl[1] && !dec_trbl[2] && !dec_trbl[3]){
            if(!(theta_p > a && theta_p < 7.0*a)){
              weighted_intensity = w;
            }
          }
          else if(dec_trbl[0] && !dec_trbl[1] && !dec_trbl[2] && !dec_trbl[3]){
            if(theta_p > a && theta_p < 3.0*a){
              weighted_intensity = w;
            }
          }
          swirl_pic.set_intensity(p.y, p.x, weighted_intensity);
        }
      }
    }
  }
  return swirl_pic;
}
