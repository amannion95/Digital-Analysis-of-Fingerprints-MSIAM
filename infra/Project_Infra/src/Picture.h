#ifndef PICTURE_H
#define PICTURE_H

#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <stdlib.h>
#include "Usefull_functions.h"
#include <algorithm>    // std::max

class Picture{
  private:

    cv::Mat picture;
    unsigned int x_length, y_length;

  public:

    Picture(const std::string& filename);
    Picture(unsigned int x_length,unsigned int y_length);
    Picture(const cv::Mat& pic);

    float get_intensity(unsigned int i, unsigned int j)const;

    void set_intensity(unsigned int i, unsigned int j,float intensity);

    void print_picture()const;

    float maximum_intensity();
    float minimum_intensity();

    Picture symmetry_transform()const;

};



























#endif
