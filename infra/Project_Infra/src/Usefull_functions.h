#ifndef USEFULL_FUNCTIONS_H
#define USEFULL_FUNCTIONS_H

#include <opencv2/opencv.hpp>



/*  ------
iitof : Takes color intensity as an int between [0,255] and converts it into a float between [0,1]
----     */
float iitof(int color);

/*  ------
iftoi : Takes color intensity a float between [0,255] and converts it into an int between [0,1]
----     */
int iftoi(float color);

/*  ------
display_mat : Takes a matrix, it's number of row, it's number of col and display it
----     */
void display_matrix(float ** matrix,int row, int col);

/* ----------
coefficient functions
---------- */
float log_coeff_isotropic(cv::Point p, cv::Point centre, double c);
//float pow_coeff_isotropic(cv::Point p, cv::Point centre, int n, float c);

// --------------------- PART 1 Max -------------------//
cv::Point rotation_ij(cv::Point a ,float rotation);
float fct_c_test(float x);
std::vector<cv::Point> segment(cv::Point a, cv::Point b);
bool compare_polar_cord(cv::Point2f a, cv::Point2f b);
bool compare_y_cord(cv::Point a, cv::Point b);

//--------------------point rotation---------------------//
cv::Point2d pt_polar_rotation(cv::Point2d orig, cv::Point2d centre, double angle);
















#endif
