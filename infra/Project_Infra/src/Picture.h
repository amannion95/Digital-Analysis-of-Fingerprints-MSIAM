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
    // add by tristan 10th jan

    Picture apply_gaussian_blur(int win_size)const;
    cv::Point get_index_maximum_intensity()const;
    cv::Point get_index_minimum_intensity()const;

    void print_pression_center(int size_win_gauss)const;
    cv::Point pressure_center_gauss();

    //neighbourhood functions
    std::vector<cv::Point> ellipse_nbh(cv::Point, unsigned int a, unsigned int b);
    void show_nbh(std::vector<cv::Point> nbh)const;

    //testing coefficient functions
    Picture log_transform_isotropic(cv::Point p, unsigned int a, unsigned int b, double c);
    Picture pow_transform_isotropic(cv::Point p, unsigned int a, unsigned int b, double c);


    //find pressure center (threshold + gaussian)
    Picture apply_threshold(float set_lim);
    std::vector<cv::Point> get_0intensity_index ();
    cv::Point get_median_center(std::vector<cv::Point> intensity_index);
    void print_median_center(int thresh);



   //return the ellipse
   Picture extract_ellipse_pic(cv::Point center, unsigned int a,unsigned int b);


   //---------------------------test--------------------------

   Picture try_apply_isotrop(cv::Point center,unsigned int a,unsigned int b);
   

};



















#endif
