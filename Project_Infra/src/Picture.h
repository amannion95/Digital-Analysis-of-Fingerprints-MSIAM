#ifndef PICTURE_H
#define PICTURE_H

#include "Useful_functions.h"
#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <cstring>


class Picture{

  private:

    cv::Mat picture;
    unsigned int x_length, y_length;

  public:
    Picture(const std::string& filename);
    Picture(unsigned int x_length,unsigned int y_length);
    Picture(const cv::Mat& pic);
    Picture();
    /**
   *\brief Function which return the intensity of the pixel (j,i).
   *
   *\param i of the coordinates (i,j) of the pixel.
   *\param j of the coordinate (i,j) of the pixel.
   *\return The intensity of the pixel (j,i).
   */
    float get_intensity(unsigned int i, unsigned int j)const;
    /**
   *\brief Function which set the intensity of the pixel (j,i).
   *
   *\param i of the coordinates (i,j) of the pixel.
   *\param j of the coordinate (i,j) of the pixel.
   *\return NONE
   */
    void set_intensity(unsigned int i, unsigned int j,float intensity);
    void operator=(Picture pic);
   /**
   *\brief Function which return the width of the image.
   *
   *\param NONE
   *\return Width of the image (unsigned int).
   */
    unsigned int get_x_len()const;
   /**
   *\brief Function which return the height of the image.
   *
   *\param NONE
   *\return Height of the image (unsigned int).
   */
    unsigned int get_y_len()const;
    /**
    *\brief Function which print the Picture.
    *
    *\param NONE
    *\return NONE
    */
    void print_picture(std::string name="Display Image")const;
    /**
   *\brief Function which return the maximum intensity of the image.
   *
   *\param NONE
   *\return Maximum intensity of the image (float in [0,1]).
   */
    float maximum_intensity()const;
    /**
   *\brief Function which return the minimum intensity of the image.
   *
   *\param NONE
   *\return Minimum intensity of the image (float in [0,1]).
   */
    float minimum_intensity()const;
    /**
   *\brief Function which perform the symetry of the image along the Y axis.
   *
   *\param NONE
   *\return The transformed Picture (type Picture).
   */
    Picture symmetry_wrt_y()const;
    /**
   *\brief Function which perform the symetry of the image along the X axis.
   *
   *\param NONE
   *\return The transformed Picture (type Picture).
   */
    Picture symmetry_wrt_x()const;
    /**
   *\brief Function which perform the symetry of the image along the top to bottom diagonal.
   *
   *\param NONE
   *\return The transformed Picture (type Picture).
   */
    Picture diagonal_symmetry_top_to_bottom()const;
    /**
   *\brief Function which perform the symetry of the image along the bottom to top diagonal.
   *
   *\param NONE
   *\return The transformed Picture (type Picture).
   */
    Picture diagonal_symmetry_bottom_to_top()const;
    /**
   *\brief Copy constructor.
   *
   *\param NONE
   *\return A copy of the Picture (type Picture).
   */
    Picture clone()const;
    /**
   *\brief Function that rescale every color such that the maximum intensity is 1 and the minimum intensity is 0.
   *
   *\param NONE
   *\return The transformed Picture (type Picture).
   */
    void rescale_color();
    /**
   *\brief Function that save the picture as a .png.
   *
   *\param name, the string that will be used as a filename (the filename has to include the extension of the image).
   *\return NONE
   */
    void SAVE_PIC(std::string name)const;
    /**
   *\brief Function that return the matrix encoding the Picture.
   *
   *\param NONE
   *\return An Eigen Matrix encoding the Picture (type: Matrix<float,Dynamic,Dynamic>).
   */
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_matrix()const;
    /**
   *\brief Function that compute the barycenter of the black pixel after a binary threshold.
   *
   *\param NONE
   *\return The barycenter of the Picture.
   */
    cv::Point center_of_pressure()const;
    // add by tristan 10th jan

    Picture apply_gaussian_blur(int win_size)const;
    cv::Point get_index_minimum_intensity()const;

    //modified 15th jan stich it on new one
    void print_pression_center_gauss_threshold();
    cv::Point pressure_center_gauss_threshold()const;

    /**
    *\brief Constructs an elliptical region on the image by checking whether or
            not the coordinates of each point satisfy a parametric equation.
    *
    *\param p Point object for the centre of the ellipse and integers for the major
            and minor axes
    *\param a minor axis of the ellipse
    *\param b major axis of the ellipse
    *
    *\return vector of Point objects corresponding to the elliptical region
    */
    std::vector<cv::Point> ellipse_nbh(cv::Point p, unsigned int a, unsigned int b)const;

    /**
    *\brief Displays the elliptical region formed by the method "ellipse_nbh" as
            a white ellipse on the image. Could be used in general as a method to
            set any given set of pixels to white.
    *
    *\param nbh vector of Point objects
    *
    *\return A copy of the image with all of the points in the input vector set to
            intensity 1.0
    */
    void show_nbh(std::vector<cv::Point> nbh)const;
    /**
   *\brief Function that perform a threshold on the image.
   *
   *\param set_lim The threshold coefficient (float in [0,1].
   *\return The transformed Picture.
   */
   Picture apply_threshold(float set_lim)const;

   //return the ellipse
   Picture extract_ellipse_pic(cv::Point& center, unsigned int a,unsigned int b)const;

   //-----------------Image Loading, Saving and Pixels Manipulation-----------//

   /**
  *\brief Function that perform an isotropic intensity transformation on the image.
  *
  *\param p the center that will be used for the isotropic function (type: Point)).
  *\return The transformed Picture.
  */
   Picture transform_isotropic(cv::Point p)const;
   /**
  *\brief Function that perform an anisotropic intensity transformation on the image using an elliptic distance given with a center, and two radius.
  *
  *\param p he center that will be used for the anisotropic function.
  *\param a the small radius of the ellipse.
  *\param b the big radius of the ellipse.
  *\return The transformed Picture
  */
   Picture transform_anisotropic(cv::Point p, unsigned int a, unsigned int b)const;
   /**
  *\brief Function that generate a random border around an ellipse
  *
  *\param center The center of the ellipse.
  *\param a the small radius of the ellipse.
  *\param b the big radius of the ellipse.
  *\return A vector containing the points of the random border generated (type: vector<Point>).
  */
   std::vector<cv::Point> weak_pressure_border(cv::Point center, unsigned int a, unsigned int b)const;
   /**
  *\brief Function that generate a random area around an ellipse.
  *
  *\param border a border around an ellipse.
  *\param center The center of the ellipse.
  *\param a the small radius of the ellipse.
  *\param b the big radius of the ellipse.
  *\return A vector containing the points of the random area generated (type: vector<Point>)
  */
   std::vector<cv::Point> weak_pressure_area(std::vector<cv::Point> border, cv::Point center, unsigned int a, unsigned int b )const;
   /**
  *\brief Function that perform an anisotropic pixel transformation on a random area around an ellipse.
  *
  *\param center The center of the ellipse.
  *\param a the small radius of the ellipse.
  *\param b the big radius of the ellipse.
  *\return The transformed Picture
  */
   Picture attenuation_weak_area(cv::Point center, unsigned int a, unsigned int b)const;



//-----------------------------ROTATION-----------------------------------------

  /**
  *\brief Checks whether a given pixel is within the frame of the image.
  *
  *\param p Point object (integer coordinates).
  *
  *\return Boolean: true if the point is within the frame, false if not.
  */
  bool isinframe(cv::Point p)const;

  /**
  *\brief Finds the barycentre - centre of pressure - of a local region of the
        image.
  *
  *\param region A list of the pixels contained in the region - Point objects.
  *
  *\return Point object representing the centre of pressure.
  */
  cv::Point local_cop(std::list<cv::Point>& region)const;

  /**
  *\brief Calls the function "local_cop" and finds the Euclidean distance from
        the returned point to another point - the centre of rotation of the
        swirl. This method was written solely to be called by "swirl_zonalmorph".
  *
  *\param subset A list of pixels to be passed to "local_cop"
  *\param centre Point object representing the centre.
  *
  *\return Double representing the Euclidean distance from the barycentre of the
         sub-region to the centre of rotation.
  *
  */
  double local_cop_distance(std::list<cv::Point> subset, cv::Point centre)const;

  /**
  *\brief Defines a square around the area on which the swirl transformation is
        being applied and divides it into grid squares. This was written solely
        to be called by the method "swirl_zonalmorph".
  *
  *\param radius The radius of the circle on which the rotation is applied (half
        the side length of the square).
  *\param subset_sidelen The number of subdivisions of the side of the square
        to make the grid squares and the centre of the
        region.
  *\param centre The centre of rotation of the swirl effect being modelled.
  *
  *\return A Matrix object (Eigen library) with lists of Point objects as entries.
  */
  Eigen::Matrix<std::list<cv::Point>, Eigen::Dynamic, Eigen::Dynamic>
  morph_subsets(int radius, int subset_sidelen, cv::Point centre)const;

  /**
  *\brief Method that performs rotation of the pixel coordinates using a
        cartesian coordinate transformation then assigns intensities by
        casting the rotated co-ordinates to integers.
  *
  *\param centre the centre of rotation as a Point object
  *\param angle The angle of rotation in degrees.
  *
  *\return A new Picture object with the rotated intensity values.
  */
  Picture cast_rotation_cart(cv::Point centre, double angle)const;

  /**
  *\brief Method that performs rotation of the pixel co-ordinates using conversion to
        polar co-ordinates then assigns intensities by casting the rotated co-ordinates
        to integers.
  *
  *\param centre the centre of rotation as a Point object
  *\param angle The angle of rotation in degrees.
  *
  *\return A new Picture object with the rotated intensity values.
  */
  Picture cast_rotation_polar(cv::Point centre, double angle)const;

  /**
  *\brief Method that performs rotation of the pixel co-ordinates using conversion to
        polar co-ordinates then assigns intensities using nearest-neighbour interpolation;
        assigning to the new pixel the intensity of the original pixel with minimal Euclidean
        distance from it.
  *
  *\param centre The centre of rotation as a Point object.
  *\param angle The angle of rotation in degrees.
  *
  *\return A new Picture object with the rotated intensity values.
  */
  Picture nn_rotation_polar(cv::Point& centre, double angle)const;

  /**
  *\brief Function that performs rotation of the pixel co-ordinates using conversion to
        polar co-ordinates then assigns intensities using bilinear interpolation; taking
        weighted sums of the surrounding pixels.
  *
  *\param centre The centre of rotation as a Point object
  *\param angle The angle of rotation in degrees.
  *
  *\return A new Picture object with the rotated intensity values.
  */

  Picture bilinear_rotation_polar(cv::Point& centre, double angle)const;

  /**
  *\brief Method that creates the swirl effect caused by rotating a finger under
        pressure in a similar way to the polar rotation functions except that in
        this case the angle is given as a function of the radial distance from
        the centre of rotation.
  *
  *\param centre The centre of the rotation, a floating-point coefficient controlling
        the "amount" of swirl to be applied.
  *\param radius The radius of the circle on which the swirl effect is to be localised.
  *
  *\return A new Picture object with the transformation applied.
  */
  Picture swirl(cv::Point centre, double twist, int radius)const;

  /**
  *\brief Pointwise bilinear interpolation: more convenient for the swirl function
        implementation.
  *
  *\param p The Point object (with double-precision floating point coordinates)
        to be interpolated.
  *
  *\return Intensity value resulting from the interpolation.
  */
  float bilinear_interpolation(cv::Point2d p)const;

  /**
  *\brief Performs morphological erosion on an elliptical region on the image.
  *
  *\param centre Point object defining the centre of the elliptical area to be
        modified
  *\param a Integer minor axis
  *\param b Integer major axis
  *
  *\return Picure object with the transformation applied.
  */
  Picture local_erosion(cv::Point centre, int a, int b)const;

  /**
  *\brief Performs morphological dilation on an elliptical region on the image.
  *
  *\param centre Point object defining the centre of the elliptical area to be modified
  *\param a Integer minor axis
  *\param b Integer major axis
  *
  *\return{Picure object with the transformation applied.}
  */
  Picture local_dilation(cv::Point centre, int a, int b)const;

  /**
  *\brief Performs the swirl transform followed by localised morphological filtering
        depending on the relative position of the barycentre of each list of
        Point objects returned by the function "morph_subsets".
  *
  *\param centre Point object defining the centre of the rotation
  *\param twist Double-precision parameter controlling the "strength" of the swirl effect, integer radius
        of the rotation.
  *\param root_num_sub Integer square root of the desired number of
        grid squares.
  *
  *\return Picture object with the transformation applied.
  */
  Picture swirl_zonalmorph(cv::Point centre, double twist, int radius, int root_num_sub)const;

  //----------------Linear Filtering----------------//


  /**
 *\brief Function that perform a discrete convolution between a Picture and a mask.
 *
 *\param mask the eigen matrix of the mask that will be convolved with the picture (type: Matrix<float,Dynamic,Dynamic>).
 *\return The result Picture of the convolution.
 */
  Picture discrete_convolution(Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> mask)const;
  /**
 *\brief Function that perform a discrete convolution between a Picture and a mask using Fast Fourrier Transform.
 *
 *\param mask the eigen matrix of the mask that will be convolved with the picture (type: Matrix<float,Dynamic,Dynamic>)
 *\return The result Picture of the convolution.
 */
  Picture ConvolutionDFT(Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> mask)const;
  /**
 *\brief Function that perform a discrete convolution between a Picture and an evolutive preset mask (given by the function Evolutive_kernel).
 *
 *\param NONE
 *\return The result Picture of the convolution.
 */
  Picture discrete_convolution_evolutive_kernel()const;
  /**
 *\brief Function that perform a discrete convolution between a Picture and an evolutive preset mask (given by the function Evolutive_kernel_no_identity).
 *
 *\param NONE
 *\return The result Picture of the convolution.
 */
  Picture discrete_convolution_evolutive_kernel_no_identity()const;


//--------------------------OPTIMIZATION PART-----------------------------------
  /**
  *\brief  Function which perform picture translation along x.
  *
  *\param  coeff translation param (must be integer)
  *\return The translated, along x axis, copied Picture
  */
   Picture translation_x(int coeff)const;
   /**
   *\brief  Function which perform picture translation along y.
   *
   *\param  coeff the translation param
   *\return The translated, along y axis, copied Picture
   */
   Picture translation_y(int coeff)const;

   /**
   *\brief  Fonction which test if two pictures are the same.
   *
   *\param  image the picture you want to be compared with.
   *\return True if the two pictures are the same, False otherwise.
   */
   bool is_same(Picture& image)const;
   /**
   *\brief  Fonction which compute the squared error between two pictures.
   *
   *\param  image the picture you want to be compared with.
   *\return The squared error between the two pictures.
   */
   float error(Picture& image)const;
   /**
   *\brief  Function performing the Bruteforce estimation of the translation
            parameter along x axis with squared error.
   *
   *\param  translated the shifted picture we are estimating shift param.
   *\return The closest integer to the translation parameter which minimize
            the error between the two pictures.
   */
   float print_loss_function_x_translation(Picture &translated)const;
   /**
   *\brief  Function performing the Bruteforce estimation of the translations
            parameters along x and y axis with squared error
            (take approximatelly 10 min).
   *
   *\param  translated the shifted picture we are estimating shift param.
   *\return The closest integer to the translation parameterS tx and ty
            which minimize the error between the two pictures.
   */
   cv::Point print_loss_function_xy_translation(Picture &translated)const;
   /**
   *\brief  Function performing a rapid estimation of the translation
            parameter along x axis by comparing the two barycenters
            with squared error.
   *
   *\param translated the shifted picture we are estimating shift param.
   *\return The translation parameterS tx and ty which minimize the error
            between the two pictures.
   */
   float loss_function_xt_by_barycenter(Picture& translated)const;
   /**
   *\brief  Function performing a rapid estimation of the translation
            parameters along x and y axis by comparing the two barycenters
            with squared error.
   *
   *\param translated the shifted picture we are estimating shift param.
   *return The integer translation parameter which minimize the error between
          the two pictures.
   */
   cv::Point loss_function_xyt_by_barycenter(Picture &translated)const;
   /**
   *\brief  Operator - is overloaded in order to give us the absolute error
            picture. If the result is a black picture then both pictures
            that were compared were the same.
   *
   *return The absolute error picture.
   */
   Picture operator-(Picture to_substract)const;
   /**
   *\brief  Fonction which compute the 2nd loss function between two pictures.
   *
   *\param  image the picture you want to be compared with.
   *\return The error between the two pictures - float.
   */
   float error_covariance_like(Picture &image)const;
   /**
   *\brief  Function performing a rapid estimation of the translation
            parameters along x and y axis by comparing the two barycenters
            with 2nd loss function.
   *
   *\param  translated the shifted picture we are estimating shift param.
   *return The integer translation parameterS which minimizes the error between
          the two pictures.
   */
   cv::Point loss_function_xyt_by_barycenter_covariance_error(Picture& translated)const;

   /**
   *\brief  Function which allow us to perform translation with floating param.
   *
   *\param  x floating shift parameter along x axis.
   *\param  y floating shift parameter along y axis.
   *return The translated picture.
   */
   Picture floating_translation(float x, float y)const;
   /**
   *\brief  Function which allow us find a more accurate translation parameter
            along x axis given the closest integer approximating shift param.
   *
   *\param  aproxim the return of one of the function loss_function_xt_by_barycenter
               or print_loss_function_x_translation.
   *\param  translated the translated picture.
   *return The floating translation parameter which minimize the error between
          the two pictures.
   */
   float find_opti_px(float aproxim,Picture &translated)const;
   /**
   *\brief  Function which allow us find more accurate translation parameters
            along x and y axis given the closest integer approximating shift param.
   *
   *\param  aproxim the returned param of one of the function loss_function_xt_by_barycenter
               or print_loss_function_x_translation.
   *\param  translated the translated picture.
   *return The floating translation parameters which minimizes the error between
          the two pictures.
   */
   cv::Point2f find_opti_px_py(cv::Point_<float> aproxim,Picture &translated)const;
   /**
   *\brief  Function which allow us to perform a fast translation with integer param
   *
   *\param x the translation param along x axis
   *\param y the translation param along y axis
   *
   *return The translated picture
   */
   Picture translation_opti_int_xy(int x,int y)const;
   /**
   *\brief  Function which approximate translations parameters along x and y axis
            using discrete Fourier transform
   *
   *\param  translated the translated picture
   *
   *return The floating translation parameters which minimize the squared error between
          the two pictures
   */
   cv::Point2f estimation_tranlsation_by_dft(Picture& translated)const;
   /**
   *\brief  Function performing the Bruteforce estimation of the rotation
            parameters with squared error
   *
   *\param  rotated the rotated picture
   *
   *return The closest integer to the rotation parameter between the two pictures
   */
   float estimation_rotation_bruteforce(Picture& rotated)const;
   /**
   *\brief  Function which allow us to accurate the rotation parameter with squared error
   *
   *\param  angle the angle approximated with estimation rotation functions
   *\param rotated the rotated picture
   *
   *return The closest float to the rotation parameter between the two pictures
   */
   float more_accurate_rotation_parameter(float angle, Picture& rotated)const;
   /**
   *\brief  Function translate the picture in order to put the barycenter of the
            picture at the center of the picture
   *
   *return New picture which have barycenter in the center
   */
   Picture put_barycenter_at_picture_center()const;
   /**
   *\brief  Function which compute the sum of the pixels intensity of a picture
   *
   *return The pixel intensity sum
   */
   float sum_intensity_picture()const;
   /**
   *\brief  Function which compute the closest interger to rotation parameter and
            the two floating points which approximate the translation parameterS
            between two pictures
   *
   *\param  Rot_txy_picture the rotated and translated picture.
   *
   *return NONE
   */
   void execution_evaluation_rtxy(Picture& Rot_txy_picture)const;


   //-------------------------------------gradient descent----------------------
   /**
   *\brief  Compute and return the partial derivative of the picture along x axis
            (finite difference scheme)
   *
   *
   *return The 2D matrix containing the intensities
   */
   float** dg_dwx()const;
   /**
   *\brief  Compute and return the partial derivative of the picture along y axis
            (finite difference scheme)
   *
   *
   *return The 2D matrix containing the intensities
   */
   float** dg_dwy()const;

   /**
   *\brief  Compute the Gradient of the squared error loss function
   *
   *\param  p_after Translation param
   *\param translated translated picture of which we are looking for translation parameters
   *
   *return The gradient of the loss function
   */
   cv::Point2f Error_partial_deriv_px_py(cv::Point2f p_after,Picture& translated)const;

   /**
   *\brief  Compute the translation param between two pictures with the gradient
            descent method. Only works if |px|+|py|<=3. Better results if both
            of them are positives
   *
   *\param  p_after translation param
   *\param  epsilon stop condition
   *\param alpha
   *\param translated the shifted picture of which we are looking for translation parameters
   *
   *return The gradient of the loss function
   */
   cv::Point2f gradient_descent(cv::Point2f p_after,float epsilon,float alpha,Picture& translated)const;


};


#endif
