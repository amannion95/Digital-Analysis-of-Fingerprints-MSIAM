#ifndef USEFUL_FUNCTIONS_H
#define USEFUL_FUNCTIONS_H

#include <opencv2/opencv.hpp>
#include <eigen3/Eigen/Dense>
#include <opencv2/core/eigen.hpp>



/**
*\brief Function that remap the intensity.
*
*\param color an integer in [0,255]. 0 corresponds to black, 255 corresponds to white.
*\return float between 0 and 1. 0 corresponds to black, 1 corresponds to white.
*/
float iitof(int color);
/**
*\brief Function that remap the intensity.
*
*\param color a float in [0,1]. 0 corresponds to black, 1 corresponds to white.
*\return int between in [0,255]. 0 corresponds to black, 255 corresponds to white.
*/
int iftoi(float color);
/**
*\brief Function that display a matrix.
*
*\param matrix a matrix of float.
*\return NULL
*/
void display_matrix(float ** matrix,int row, int col);
/**
*\brief Function that return the Points given by the segment [a,b].
*
*\param a starting Point of the segment [a,b].
*\param b  end Point of the segment [a,b].
*\return a float between 0 and 1. 0 corresponds to black, 1 corresponds to white.
*/
std::vector<cv::Point> segment(cv::Point a, cv::Point b);
/**
*\brief Order relation between two Points a and b (polar coordinate). Point are ordered by angle, then by norm.
*
*\param a first Point to compare.
*\param b second Point to compare.
*\return a boolean. True if a<b, false if not.
*/
bool compare_polar_cord(cv::Point2f a, cv::Point2f b);
/**
*\brief Order relation between two Points a and b (carthesian coordinate). Point are ordered by y coordinate, then x coordinate.
*
*\param a first Point to compare.
*\param b second Point to compare.
*\return a boolean. True if a<b, false if not.
*/
bool compare_y_cord(cv::Point a, cv::Point b);
/**
*\brief perform the convolution (using FFT) of A and B and return it in C.
*
*\param A first matrix to convolve.
*\param B second matrix to convolve.
*\param C a matrix that will store the result of the convolution.
*\return NULL
*/
void convolveDFT(cv::InputArray A, cv::InputArray B, cv::OutputArray C);
/**
*\brief Function that compute an evolutive kernel matrix which is the weighted mean of an evolutive gaussain blurr and the identity kernel.
*
*\param t is the weight parameter.
*\param size is the size of the gaussian kernel tu use.
*\return a kernel matrix
*/
Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Evolutive_kernel(float t, int size);
/**
*\brief Function that compute an evolutive gaussian blurr
*
*\param t is the weight parameter.
*\param size is the size of the gaussian kernel tu use.
*\return a kernel matrix
*/
Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Evolutive_kernel_no_identity(float t, int size);

/**
brief Rotates a single point through a given angle using conversion to polar coordinates
*
*\param orig Point to be rotated
*\param centre Centre of rotation
*\param angle Angle through which to rotate the point
*
*\return Point object containing the coordinates the original point would have in a
        rotated image.
*/
cv::Point2d pt_polar_rotation(cv::Point2d orig, cv::Point2d centre, double angle);




















#endif
