#include <iostream>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <math.h>
#include <algorithm>
#include <iterator>

#include "Usefull_functions.h"
#include "Picture.h"

using namespace std;
using namespace cv;
using namespace Eigen;

Point Picture::local_cop(list<Point> region){
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

Matrix<list<Point>, Dynamic, Dynamic> Picture::morph_subsets(int radius, int subset_sidelen, Point centre){
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

double Picture::local_cop_distance(list<Point> subset, Point centre){
  Point cop = local_cop(subset);
  double dist = norm(cop - centre);
  return dist;
}

Picture Picture::swirl_morph(Point centre, double twist, int radius, int subset_sidelen){
  Picture swirl_pic(x_length, y_length);
  //cout << "subset 1: " << subsets[0] << endl;
  for(int i=0; i<x_length; i++){
    for(int j=0; j<y_length; j++){
      double dist = norm(centre - Point(i,j));
      if(dist < radius){
        double amt = 1 - dist / (double)radius;
        double angle = 2 * twist * M_PI * amt;
        //cout << angle << endl;
        Point2d p = pt_polar_rotation(Point2d((double)i, (double)j), centre, -1*angle);
        float intensity = this->bilinear_interpolation(p);
        swirl_pic.set_intensity(j, i, intensity);
        /*if(intensity != 1.0){
          cout << intensity << endl;
        }*/
      }
      else{
        swirl_pic.set_intensity(j, i, get_intensity(j, i));
      }
    }
  }
  cout << "swirled!" << endl;
  Matrix<list<Point>, Dynamic, Dynamic> subsets = this->morph_subsets(radius, subset_sidelen, centre);
  int c = 0;
  int d = 0;
  int e = 0;
  int f = 0;
  for(size_t j=0; j<subsets.cols(); ++j){
    for(size_t i=0; i<subsets.rows(); ++i){
      //cout << "matrix loopin" << endl;
      list<Point> subset = subsets(i,j);
      Point topl = subset.front();
      Point botr = subset.back();
      if(isinframe(topl) && isinframe(botr)){
        //cout << "into morph" << endl;
        float d1 = this->local_cop_distance(subset, centre);
        float d2 = swirl_pic.local_cop_distance(subset, centre);
        //cout << "got cops: " << d1 << " and " << d2 << endl;
        if(d1 > d2){
          e++;
          //local erosion on subset
          swirl_pic.bounded_erosion(subset);
          /*for(Point& p : subset){
            swirl_pic.set_intensity(p.y, p.x, 0.0);
          }*/
        }
        else if(d1 < d2){
          d++;
          //local dilation on subset
          swirl_pic.bounded_dilation(subset);
          /*for(Point& p : subset){
            swirl_pic.set_intensity(p.y, p.x, 1.0);
          }*/
        }
        else{c++;}
      }
      else{f++;}
    }
  }
  cout << e << " erosions\n" << d << " dilations\n" << c << " neither\n" << f << " oob" << endl;
  return swirl_pic;
}

void Picture::bounded_erosion(list<Point> points){
  /*for(Point& p : points){
    if(p.x < 0 && p.x >= x_length && p.y < 0 && p.y >= y_length){
      cout << "check 1 " << p <<  endl;
      points.erase(remove(points.begin(), points.end(), p), points.end());
      cout << "check 2 " << p << endl;
    }
  }
  cout << "frame sorted" << endl;*/
  list<int> x_list, y_list;
  for(Point& p : points){
    x_list.emplace_back(p.x);
    y_list.emplace_back(p.y);
  }
  int xmax = *max_element(x_list.begin(), x_list.end());
  int ymax = *max_element(y_list.begin(), y_list.end());
  int xmin = *min_element(x_list.begin(), x_list.end());
  int ymin = *min_element(y_list.begin(), y_list.end());
  int halfsize = (int)round((xmax-xmin)/2);
  cout << "e:half subset size: " << halfsize << endl;
  int in = 0;
  int mid = 0;
  int out = 0;
  Point sub_cen(xmin+halfsize, ymin+halfsize);
  cout << "centre: " << sub_cen << endl;
  Picture temp_pic = clone();
  for(Point& p : points){
    double d = norm(p-sub_cen);
    vector<float> elem;
    cout << "distance: " << d <<  endl;
    if(d<=halfsize*0.5){
      cout << "inner loop" << endl;
      in++;
      for(int i=-1; i<=1; i++){
        //for(int j=-1; j<=1; j++){
        elem.push_back(temp_pic.get_intensity(p.y, p.x+i));
      }
    }
    else if(d>halfsize*0.5 && d<=halfsize*0.75){
      cout << "middle loop" << endl;
      mid++;
      //double r = norm(p-cop);
      for(int i=-1; i<=1; i++){
        elem.push_back(temp_pic.get_intensity(p.y, p.x+i));
      }
      for(int j=-1; j<=1; j++){
        elem.push_back(temp_pic.get_intensity(p.y+j, p.x));
      }
    }
    else{
      out++;
      cout << "outer loop" <<  endl;
    }
    float min = *min_element(elem.begin(), elem.end());
    set_intensity(p.y, p.x, min);
    elem.clear();
  }
  cout << "inner:" << in << ", middle:" << mid << ", outer:" << out <<
  ", total: " << in+mid+out << endl;
}

void Picture::bounded_dilation(list<Point> points){
  /*for(Point& p : points){
    if(p.x < 0 && p.x >= x_length && p.y < 0 && p.y >= y_length){
      cout << "check 1 " << p <<  endl;
      points.erase(remove(points.begin(), points.end(), p), points.end());
      cout << "check 2 " << p << endl;
    }
  }
  cout << "frame sorted" << endl;*/
  list<int> x_list, y_list;
  for(Point& p : points){
    x_list.emplace_back(p.x);
    y_list.emplace_back(p.y);
  }
  int xmax = *max_element(x_list.begin(), x_list.end());
  int ymax = *max_element(y_list.begin(), y_list.end());
  int xmin = *min_element(x_list.begin(), x_list.end());
  int ymin = *min_element(y_list.begin(), y_list.end());
  int halfsize = (int)round((xmax-xmin)/2);
  Point sub_cen(xmin+halfsize, ymin+halfsize);
  //cout << "d:half subset size: " << halfsize << endl;
  int in = 0;
  int mid = 0;
  int out = 0;
  Picture temp_pic = clone();
  for(Point& p : points){
    double d = norm(p - sub_cen);
    vector<float> elem;
    if(d<=halfsize*0.5){
      in++;
      for(int i=-1; i<=1; i++){
        for(int j=-1; j<=1; j++){
          elem.push_back(temp_pic.get_intensity(p.y, p.x+i));
        }
      }
    }
    else if(d>halfsize*0.5 && d<=halfsize*0.75){
      mid++;
      for(int i=-1; i<=1; i++){
        elem.push_back(temp_pic.get_intensity(p.y, p.x+i));
      }
      for(int j=-1; j<=1; j++){
        elem.push_back(temp_pic.get_intensity(p.y+j, p.x));
      }
    }
    else{out++;}
    float max = *max_element(elem.begin(), elem.end());
    set_intensity(p.y, p.x, max);
    elem.clear();
  }
  cout << "inner:" << in << ", middle:" << mid << ", outer:" << out <<
  ", total:" << in+mid+out << endl;
}

Picture Picture::swirl_erode(Point centre, double twist, int radius){
  Picture swirl_pic(x_length, y_length);
  for(int i=0; i<x_length; i++){
    for(int j=0; j<y_length; j++){
      double dist = norm(centre - Point(i,j));
      if(dist < radius){
        double amt = 1 - dist / (double)radius;
        double angle = 2 * twist * M_PI * amt;
        Point2d p = pt_polar_rotation(Point2d((double)i, (double)j), centre, -1*angle);
        float intensity = this->bilinear_interpolation(p);
        swirl_pic.set_intensity(j, i, intensity);
      }
      else{
        swirl_pic.set_intensity(j, i, get_intensity(j, i));
      }
    }
  }
  list<Point> erosionbox;
  int zero_x = centre.x - radius;
  int zero_y = centre.y - radius;
  cout << "zero: " << Point(zero_x, zero_y) << endl;
  for(unsigned int i=zero_x; i<=zero_x+2*radius; i++){
    for(unsigned int j=zero_y; j<=zero_y+2*radius; j++){
      erosionbox.emplace_back(Point(i,j));
    }
  }
  swirl_pic.bounded_erosion(erosionbox);
  return swirl_pic;
}

Picture Picture::swirl_dilate(Point centre, double twist, int radius){
  Picture swirl_pic(x_length, y_length);
  for(int i=0; i<x_length; i++){
    for(int j=0; j<y_length; j++){
      double dist = norm(centre - Point(i,j));
      if(dist < radius){
        double amt = 1 - dist / (double)radius;
        double angle = 2 * twist * M_PI * amt;
        Point2d p = pt_polar_rotation(Point2d((double)i, (double)j), centre, -1*angle);
        float intensity = this->bilinear_interpolation(p);
        swirl_pic.set_intensity(j, i, intensity);
      }
      else{
        swirl_pic.set_intensity(j, i, get_intensity(j, i));
      }
    }
  }
  list<Point> dilationbox;
  int zero_x = centre.x - radius;
  int zero_y = centre.y - radius;
  for(unsigned int i=zero_x; i<=zero_x+2*radius; i++){
    for(unsigned int j=zero_y; j<=zero_y+2*radius; j++){
      dilationbox.emplace_back(Point(i,j));
    }
  }
  swirl_pic.bounded_dilation(dilationbox);
  return swirl_pic;
}

/*Picture Picture::swirl_morph_zones(Point centre, double twist, int radius, int subset_sidelen){
  Picture swirl_pic(x_length, y_length);
  for(int i=0; i<x_length; i++){
    for( int j=0; j<y_length; j++){
      double dist = norm(centre - Point(i,j));
      if(dist < radius){
        double amt = 1 - dist / (double)radius;
        double angle = 2 * twist * M_PI * amt;
        Point2d p = pt_polar_rotation(Point2d((double)i, (double)j), centre, -1*angle);
        float intensity = this->bilinear_interpolation(p);
        swirl_pic.set_intensity(j, i, intensity);
      }
      else{
        swirl_pic.set_intensity(j, i, get_intensity(j, i));
      }
    }
  }
  Matrix<list<Point>, Dynamic, Dynamic> subsets = this->morph_subsets(radius, subset_sidelen, centre);
  MatrixXi morph_indic(subset_sidelen, subset_sidelen);
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
      else{}
    }
  }
  for(size_t j=0; j<subsets.cols(); j++){
    for(size_t i=0; i<subsets.rows(); i++){
      list<int> x_list, y_list;
      for(Point& p : subsets(i,j)){
        x_list.emplace_back(p.x);
        y_list.emplace_back(p.y);
      }
      int xmax = *max_element(x_list.begin(), x_list.end());
      int ymax = *max_element(y_list.begin(), y_list.end());
      int xmin = *min_element(x_list.begin(), x_list.end());
      int ymin = *min_element(y_list.begin(), y_list.end());
      int halfsize = (int)round((xmax-xmin)/2);
      Point sub_cen(xmin+halfsize, ymin+halfsize);
      //int halfarea = 2 * halfsize * halfsize;
      Matrix<Point, 4, Dynamic> twoquads; //clockwise!!
      /*Matrix<Point, 3, Dynamic> top_edge, right_edge, bottom_edge, left_edge;
      Matrix<Point, 2, Dynamic> topleft, topright, bottomright, bottomleft;
      for(unsigned int i=xmin; i<=xmax; i++){
        for(unsigned int j=sub_cen.y; j<=ymax; j++){
          twoquads(0, i*(ymax-sub_cen.y)+j) = Point(i,j);
        }
      }
      for(unsigned int i=sub_cen.x; i<=xmax; i++){
        for(unsigned int j=ymin; j<=ymax; j++){
          twoquads(1, i*(ymax-ymin)+j) = Point(i,j);
        }
      }
      for(unsigned int i=xmin; i<=xmax; i++){
        for(unsigned int j=ymin; j<=sub_cen.y; j++){
          twoquads(2, i*(sub_cen.y-ymin)+j) = Point(i,j);
        }
      }
      for(unsigned int i=xmin; i<=sub_cen.x; i++){
        for(unsigned int j=ymin; j<=ymax; j++){
          twoquads(3, i*(ymax-ymin)+j) = Point(i,j);
        }
      }
      switch(morph_indic(i,j)){
        case 1 :
          swirl_pic.check_nbh_erosion(morph_indic, i, j, twoquads, sub_cen);
          break;
        case -1 :
          swirl_pic.check_nbh_dilation(morph_indic, i, j, twoquads, sub_cen);
          break;
        default :
          break;
      }
    }
  }
}*/

/*void Picture::check_nbh_erosion(MatrixXi morph_indic, unsigned int i, unsigned int j, Matrix<Point, 4, Dynamic> twoquads){
  //list<Point> tophalf, righthalf, bottomhalf, lefthalf;
  MatrixXi padded(morph_indic.rows()+2, morph_indic.cols()+2);
  for(size_t l=0; l<=padded.cols(); l++){
    padded(0,l) = 0;
    padded(padded.rows(),l) = 0;
  }
  for(size_t k=0; k<=padded.rows(); k++){
    padded(k,0) = 0;
    padded(k,padded.cols()) = 0;
  }
  for(size_t l=0; l<=morph_indic.cols(); l++){
    for(size_t k=0; k<=morph_indic.rows(); k++){
      padded(k+1,l+1) = morph_indic(k,l);
    }
  }
  /*for(size_t l=0; l<=twoquads.cols(); l++){
    tophalf.emplace_back(twoquads(0,l));
    righthalf.emplace_back(twoquads(1,l));
    bottomhalf.emplace_back(twoquads(2,l));
    lefthalf.emplace_back(twoquads(3,l));
  }
  int neighbours[4] = {padded(i+1,j), padded(i-1,j), padded(i,j-1), padded(i,j+1)};
  for(unsigned int i=0; i<4; i++){
    list<Point> h;
    for(size_t l=0; l<=twoquads.cols(); l++){
      h.emplace_back(twoquads(i,l));
    }
    if(neighbours[i]==morph_indic(i,j)){
      //constant erosion on h
    }
    else{
      //decreasing erosion on h
    }
    h.clear();
  }
}

void Picture::check_nbh_dilation(MatrixXi morph_indic, unsigned int i, unsigned int j, Matrix<Point, 4, Dynamic> twoquads){
  MatrixXi padded(morph_indic.rows()+2, morph_indic.cols()+2);
  for(size_t l=0; l<=padded.cols(); l++){
    padded(0,l) = 0;
    padded(padded.rows(),l) = 0;
  }
  for(size_t k=0; k<=padded.rows(); k++){
    padded(k,0) = 0;
    padded(k,padded.cols()) = 0;
  }
  for(size_t l=0; l<=morph_indic.cols(); l++){
    for(size_t k=0; k<=morph_indic.rows(); k++){
      padded(k+1,l+1) = morph_indic(k,l);
    }
  }
  int neighbours[4] = {padded(i+1,j), padded(i-1,j), padded(i,j-1), padded(i,j+1)};
  for(unsigned int i=0; i<4; i++){
    list<Point> h;
    for(size_t l=0; l<=twoquads.cols(); l++){
      h.emplace_back(twoquads(i,l));
    }
    if(neighbours[i]==morph_indic(i,j)){
      //constant dilation on h
    }
    else{
      //decreasing dilation on h
    }
    h.clear();
  }
}*/


Picture Picture::swirl_morph_zones(Point centre, double twist, int radius, int subset_sidelen){
  Picture swirl_pic(x_length, y_length);
  for(int i=0; i<x_length; i++){
    for( int j=0; j<y_length; j++){
      double dist = norm(centre - Point(i,j));
      if(dist < radius){
        double amt = 1 - dist / (double)radius;
        double angle = 2 * twist * M_PI * amt;
        Point2d p = pt_polar_rotation(Point2d((double)i, (double)j), centre, -1*angle);
        float intensity = this->bilinear_interpolation(p);
        swirl_pic.set_intensity(j, i, intensity);
      }
      else{
        swirl_pic.set_intensity(j, i, get_intensity(j, i));
      }
    }
  }
  Matrix<list<Point>, Dynamic, Dynamic> subsets = this->morph_subsets(radius, subset_sidelen, centre);
  //cout << "matrix width: " << subsets.cols() << endl;
  MatrixXi morph_indic(subset_sidelen, subset_sidelen);
  for(size_t j=0; j<subsets.cols(); j++){
    for(size_t i=0; i<subsets.rows(); i++){
      list<Point> subset = subsets(i,j);
      Point topl = subset.front();
      Point botr = subset.back();
      if(isinframe(topl) && isinframe(botr)){
        float d1 = this->local_cop_distance(subset, centre);
        float d2 = swirl_pic.local_cop_distance(subset, centre);
        if(d1 > d2){
          //cout << "dilation_indic_" << Point(i,j) << endl;
          morph_indic(i,j) = 1;
        }
        else if(d1 < d2){
          //cout << "erosion_indic_" << Point(i,j) << endl;
          morph_indic(i,j) = -1;
        }
        else{
          //cout << "nada_indic_" << Point(i,j) << endl;
          morph_indic(i,j) = 0;
        }
      }
      else{
        //cout << "oob_indic_" << Point(i,j) << endl;
        morph_indic(i,j) = 0;
      }
    }
  }
  cout << "morphology indicator: \n" << morph_indic << endl;
  for(size_t j=0; j<subsets.cols(); j++){
    for(size_t i=0; i<subsets.rows(); i++){
      list<int> x_list, y_list;
      for(Point& p : subsets(i,j)){
        x_list.emplace_back(p.x);
        y_list.emplace_back(p.y);
      }
      int xmax = *max_element(x_list.begin(), x_list.end());
      int ymax = *max_element(y_list.begin(), y_list.end());
      int xmin = *min_element(x_list.begin(), x_list.end());
      int ymin = *min_element(y_list.begin(), y_list.end());
      int halfsize = (int)round((xmax-xmin)/2);
      Point sub_cen(xmin+halfsize, ymin+halfsize);
      //cout << "subset centre: " << sub_cen << endl;
      swirl_pic.neighbourcheck_morph(subsets(i,j), morph_indic, i, j, sub_cen);
    }
  }
  return swirl_pic;
}

void Picture::neighbourcheck_morph(list<Point> subset, MatrixXi morph_indic, unsigned int i, unsigned int j, Point centre){
  MatrixXi padded(morph_indic.rows()+2, morph_indic.cols()+2);
  for(size_t l=0; l<padded.cols(); l++){
    padded(0,l) = 0;
    padded(padded.rows()-1,l) = 0;
  }
  for(size_t k=0; k<padded.rows(); k++){
    padded(k,0) = 0;
    padded(k,padded.cols()-1) = 0;
  }
  for(size_t l=1; l<padded.cols()-1; l++){
    for(size_t k=1; k<padded.rows()-1; k++){
      padded(k,l) = morph_indic(k-1,l-1);
    }
  }
  //cout << "padded matrix: " << padded << endl;
  //right, left, bottom, top
  int x_neighbours[2] = {padded(i,j+2), padded(i,j)};
  int y_neighbours[2] = {padded(i+2,j), padded(i,j)};
  bool x_desc[2];
  bool y_desc[2];
  for(unsigned int i=0; i<2; i++){
    if(x_neighbours[i]!=morph_indic(i,j)){
      x_desc[i] = 1;
    }
    else{
      x_desc[i] = 0;
    }
  }
  for(unsigned int i=0; i<2; i++){
    if(y_neighbours[i]!=morph_indic(i,j)){
      y_desc[i] = 1;
    }
    else{
      y_desc[i] = 0;
    }
  }
  if(morph_indic(i,j) != 0){
    zonal_morphing(x_desc, y_desc, subset, centre, morph_indic(i,j));
  }
}

/*void Picture::zonal_erosion(bool* x_desc, bool* y_desc, list<Point> subset,
  Point sub_cen){
  int x_cen = sub_cen.x;
  int y_cen = sub_cen.y;
  Picture temp_pic = clone();
  for(Point& p : subset){
    Vector3f x_coefs, y_coefs;
    int x = p.x;
    int y = p.y;
    int dx = x-x_cen;
    int dy = y-y_cen;
    if(x_desc[0] && x_desc[1]){
      for(int i=-1; i<2; i++){
        int d = (dx+i)/15;
        x_coefs(i+1) = 1-exp_morph_coef(d);
      }
    }
    else if(x_desc[0] && !x_desc[1]){//decreasing to the right, constant to the left
      if(dx>=0){
        for(int i=-1; i<2; i++){
          int d = (dx+i)/15;
          x_coefs(i+1) = 1-exp_morph_coef(d);
        }
      }
      else{
        x_coefs = {1,1,1};
      }
    }
    else if(!x_desc[0] && x_desc[1]){//constant to the right, decreasing to the left
      if(dx<=0){
        for(int i=-1; i<2; i++){
          int d = (dx+i)/15;
          x_coefs(i+1) = 1-exp_morph_coef(d);
        }
      }
      else{
        x_coefs = {1,1,1};
      }
    }
    else{
      x_coefs = {1,1,1};
    }
    if(y_desc[0] && y_desc[1]){
      for(int i=-1; i<2; i++){
        int d = (dx+i)/15;
        y_coefs(i+1) = 1-exp_morph_coef(d);
      }
    }
    else if(y_desc[0] && !y_desc[1]){//decreasing to the right, constant to the left
      if(dy>=0){
        for(int i=-1; i<2; i++){
          int d = (dx+i)/15;
          y_coefs(i+1) = 1-exp_morph_coef(d);
        }
      }
      else{
        y_coefs = {1,1,1};
      }
    }
    else if(!y_desc[0] && y_desc[1]){//constant to the right, decreasing to the left
      if(dy<=0){
        for(int i=-1; i<2; i++){
          int d = (dx+i)/15;
          y_coefs(i+1) = 1-exp_morph_coef(d);
        }
      }
      else{
        y_coefs = {1,1,1};
      }
    }
    else{
      x_coefs = {1,1,1};
    }
    Matrix3f elem = x_coefs * y_coefs.transpose();
    vector<float> intensities;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        intensities.push_back(abs(temp_pic.get_intensity(y-1+j,x-1+i)*elem(i,j)));
      }
    }
    float min = *min_element(intensities.begin(), intensities.end());
    set_intensity(y, x, min);
  }
}*/

void Picture::zonal_morphing(bool* x_desc, bool* y_desc, list<Point> subset,
  Point sub_cen, int indicator){
  int x_cen = sub_cen.x;
  int y_cen = sub_cen.y;
  int max_dist_x = abs(subset.back().x - x_cen);
  int max_dist_y = abs(subset.back().y - y_cen);
  Picture temp_pic = clone();
  for(Point& p : subset){
    //Matrix<float, 5, 1> x_coefs, y_coefs;
    Vector3f x_coefs, y_coefs;
    int x = p.x;
    int y = p.y;
    int dx = x-x_cen;
    int dy = y-y_cen;
    if(x_desc[0] && x_desc[1]){
      for(int i=-1; i<2; i++){
        int d = (dx+i)/max_dist_x;
        x_coefs(i+1) = se_coef_pow(d);
      }
    }
    else if(x_desc[0] && !x_desc[1]){//decreasing to the right, constant to the left
      if(dx>=0){
        for(int i=-1; i<2; i++){
          int d = (dx+i)/max_dist_x;
          x_coefs(i+1) = se_coef_pow(d);
        }
      }
      else{
        x_coefs = {1,1,1};
      }
    }
    else if(!x_desc[0] && x_desc[1]){//constant to the right, decreasing to the left
      if(dx<=0){
        for(int i=-1; i<2; i++){
          int d = (dx+i)/max_dist_x;
          x_coefs(i+1) = se_coef_pow(d);
        }
      }
      else{
        x_coefs = {1,1,1};
      }
    }
    else{
      x_coefs = {1,1,1};
    }
    if(y_desc[0] && y_desc[1]){
      for(int i=-1; i<2; i++){
        int d = (dy+i)/max_dist_y;
        y_coefs(i+1) = se_coef_pow(d);
      }
    }
    else if(y_desc[0] && !y_desc[1]){//decreasing to the right, constant to the left
      if(dy>=0){
        for(int i=-1; i<2; i++){
          int d = (dy+i)/max_dist_y;
          y_coefs(i+1) = se_coef_pow(d);
        }
      }
      else{
        y_coefs = {1,1,1};
      }
    }
    else if(!y_desc[0] && y_desc[1]){//constant to the right, decreasing to the left
      if(dy<=0){
        for(int i=-1; i<2; i++){
          int d = (dy+i)/max_dist_y;
          y_coefs(i+1) = se_coef_pow(d);
        }
      }
      else{
        y_coefs = {1,1,1};
      }
    }
    else{
      x_coefs = {1,1,1};
    }
    //cout << "x coeffs: " << x_coefs.transpose() << "\ny coeffs: " << y_coefs.transpose() << endl;
    Matrix3f elem = x_coefs * y_coefs.transpose();
    //cout << "product: \n" << elem << endl;
    float energy = elem.sum();
    if(energy != 0){elem /= energy;}
    //cout << "normalised: \n" << elem << endl;
    //MatrixXf elem(5,5);
    //elem = x_coefs * y_coefs.transpose();
    //cout << "Structuring element\n" << elem << endl;
    Matrix3f intensities;
    float chosen_intensity;
    if(indicator==1){
      for(int i=0; i<=2; i++){
        for(int j=0; j<=2; j++){
          float dilated_intensity = abs(temp_pic.get_intensity(y-1+j,x-1+i)+elem(i,j));
          float real_intensity;
          if(dilated_intensity > 1.0){
            real_intensity = 1.0;
          }
          else{real_intensity = dilated_intensity;}
          intensities(i,j) = real_intensity;
        }
      }
      /*float energy = intensities.sum();
      intensities /= energy;*/
      chosen_intensity = intensities.maxCoeff();
    }
    else{
      for(int i=0; i<=2; i++){
        for(int j=0; j<=2; j++){
          intensities(i,j) = (abs(temp_pic.get_intensity(y-1+j,x-1+i)-elem(i,j)));
        }
      }
      /*float energy = intensities.sum();
      intensities /= energy;*/
      chosen_intensity = intensities.minCoeff();
    }
    /*if(chosen_intensity > 1){
      cerr << "x coeffs: " << x_coefs.transpose() << "\ny_coeffs: " << y_coefs.transpose()
        << "\nnormalised product: \n" << elem << "\nintensities: \n" << intensities << endl;
    }*/
    set_intensity(y, x, chosen_intensity);
  }
}
