#ifndef PICTURE_H
#define PICTURE_H

#include "Usefull_functions.h"
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <cstdlib> // absolute value
#include <algorithm>

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

    void SAVE_PIC(std::string name);
    float** get_matrix();

    cv::Point center_of_pressure();
    // add by tristan 10th jan

    Picture apply_gaussian_blur(int win_size)const;
    cv::Point get_index_minimum_intensity()const;

    //modified 15th jan stich it on new one
    void print_pression_center_gauss_threshold();
    cv::Point pressure_center_gauss_threshold();

    //neighbourhood functions
    std::vector<cv::Point> ellipse_nbh(cv::Point, unsigned int a, unsigned int b);
    void show_nbh(std::vector<cv::Point> nbh)const;

    //testing coefficient functions
    Picture log_transform_isotropic(cv::Point p, unsigned int a, unsigned int b, double c);
    Picture pow_transform_isotropic(cv::Point p, unsigned int a, unsigned int b, double c);


    //find pressure center (threshold + gaussian)
    Picture apply_threshold(float set_lim);
    std::vector<cv::Point> get_0intensity_index ();




   //return the ellipse
   Picture extract_ellipse_pic(cv::Point center, unsigned int a,unsigned int b);


   //---------------------------test--------------------------
   //return ellipse with the right color
   Picture apply_anisotrope(cv::Point center,unsigned int a,unsigned int b);




   Picture without_noise();

   //good filter for treating ellipse after function c
   Picture accentuation_diff(int winsize );





//--------------------------------OPTIMIZATION PART-----------------------------
   Picture translation_x(int coeff);
   Picture translation_y(int coeff);
   bool is_same(Picture image);
   float error(Picture& image);
   float print_loss_function_x_translation(Picture &translated);

   cv::Point print_loss_function_xy_translation(Picture &translated);

   float loss_function_xt_by_barycenter(Picture& translated);
   cv::Point loss_function_xyt_by_barycenter(Picture &translated);


   Picture operator-(Picture to_substract);

   float error_covariance_like(Picture &image);
   cv::Point loss_function_xyt_by_barycenter_covariance_error(Picture translated);
   Picture floating_translation(float x, float y);
   float sum_error(Picture to_be_compared);
   float find_opti_px(float aproxim,Picture &translated);
   cv::Point2f find_opti_px_py(cv::Point_<float> aproxim,Picture &translated);
   Picture translation_opti_int_xy(int x,int y);
};



















#endif
