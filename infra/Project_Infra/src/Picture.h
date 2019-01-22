#ifndef PICTURE_H
#define PICTURE_H

#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <stdlib.h>
#include "Usefull_functions.h"
#include <math.h>
#include <vector>
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
    void SAVE_PIC(std::string name);
    void rescale_color();
    float** apply_c(float** c,int x_len,int y_len);
    float** get_matrix();
    std::vector<cv::Point> ellipse_nbh(cv::Point p, unsigned int a, unsigned int b);
    Picture operator-(Picture to_substract);
    void show_nbh(std::vector<cv::Point> nbh)const;
    Picture log_transform_isotropic(cv::Point p, unsigned int a, unsigned int b, double coef);
    Picture pow_transform_isotropic(cv::Point p, unsigned int a, unsigned int b, double coef);


    cv::Point center_of_pressure();
    //Rotation.cpp methods:
    std::vector<cv::Point> original_coordt(float o);
    std::vector<cv::Point_<float>> list_send_point(float o);
    Picture bilinear_interpolation_cart( double o);
    Picture cast_rotation_cart(cv::Point centre,float o);
    std::vector<cv::Point> neighb_original_coordt(float o,float a, float b, cv::Point centre);
    std::vector<cv::Point_<float>> neighb_ellipse_point(float o,float a, float b, cv::Point centre);
    Picture bilinear_interpolation_cart_nghb( double o, cv::Point centre, float a, float b);
    std::vector<cv::Point_<float>> swirl_point(float o,float a, float b, cv::Point centre);
    std::vector<cv::Point> original_swirl_point(float o,float a, float b, cv::Point centre);
    Picture bilinear_interpolation_cart_swirl( double o, cv::Point centre, float a, float b);
    Picture floating_translation(float o,float tx, float ty);
    Picture erosion_bounded(std::vector<cv::Point> v);
    Picture dilation_bounded(std::vector<cv::Point> v);
    std::vector<cv::Point> donut(cv::Point centre, int a, int b);
};


































#endif
