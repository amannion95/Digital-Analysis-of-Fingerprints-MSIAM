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
    Picture();
    float get_intensity(unsigned int i, unsigned int j)const;
    void set_intensity(unsigned int i, unsigned int j,float intensity);
    void operator=(Picture pic);
    unsigned int get_x_len();
    unsigned int get_y_len();
    void print_picture()const;
    float maximum_intensity()const;
    float minimum_intensity()const;
    Picture symmetry_wrt_y()const;
    Picture symmetry_wrt_x()const;
    Picture diagonal_symmetry_top_to_bottom()const;
    Picture diagonal_symmetry_bottom_to_top()const;
    Picture clone()const;
    void rescale_color();
    float** get_matrix();
    cv::Point center_of_pressure();
};




























#endif
