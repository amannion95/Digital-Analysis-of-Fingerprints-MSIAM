#ifndef PICTURE_H
#define PICTURE_H

#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <stdlib.h>

class Picture{
  private:
    cv::Mat picture;
    unsigned int x_length, y_length;

  public:
    Picture(const std::string& filename);
    Picture(const cv::Mat& pic);
    
    float get_intensity(unsigned int i, unsigned int j)const;

    void set_intensity(unsigned int i, unsigned int j,float intensity);

};



























#endif
