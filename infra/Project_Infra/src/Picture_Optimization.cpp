#include "Picture.h"
#include "Usefull_functions.h"

#include <time.h>
#include <fstream>

using namespace cv;
using namespace std;


Picture Picture::translation_x(int coeff){
  Picture translated=clone();
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      if ((i-coeff>=0)&&(i-coeff<x_length)){
        translated.set_intensity(j,i,get_intensity(j,i-coeff));
      }
      else{
        translated.set_intensity(j,i,1);
      }
    }
  }
  return translated;
}

Picture Picture::translation_y(int coeff){
  Picture translated=clone();
  for(int j=0;j<y_length;j++){
    for(int i=0;i<x_length;i++){
      if ((j-coeff>=0)&&(j-coeff<y_length)){
        translated.set_intensity(j,i,get_intensity(j-coeff,i));
      }
      else{
        translated.set_intensity(j,i,1);
      }
    }
  }
  return translated;
}

//fast method to do translations. Be carefull it only works for integer values.
//it copy the rectangle which will be remaining after translation then paste it
//in the right place on a white picture
Picture Picture::translation_opti_int_xy(int x,int y){
  int offset_x=max(0,-x);
  int offset_y=max(0,-y);

  Rect roi;

  roi.x = offset_x;
  roi.y = offset_y;
  roi.width = x_length-abs(x);
  roi.height = y_length-abs(y);

  Mat crop =picture(roi);
  Mat white(y_length,x_length,CV_8UC1,255);

  crop.copyTo(white(Rect(max(0,x),max(0,y),crop.cols,crop.rows)));
  Picture translated(white);
  return translated;
}

//Allow us to make a translations with floating number as parameter.
Picture Picture::floating_translation(float tx, float ty){
  Mat image2;
  Mat M=(Mat_<float>(2,3)<<1,0,tx,0,1,ty);
  warpAffine( picture, image2, M, picture.size());
  return Picture(image2);
}

//True if both pictures are the same. False otherwise
bool Picture::is_same(Picture image){
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      if (get_intensity(j,i)!=image.get_intensity(j,i)){
        return false;
      }
    }
  }
  return true;
}

//if an intensity is negative we take its abs value
//Allow us to print absolute error image.
Picture Picture::operator-(Picture to_substract){
  float ab;
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      ab=abs(get_intensity(j,i)-to_substract.get_intensity(j,i));
      to_substract.set_intensity(j,i,ab);
    }
  }
  return to_substract;
}

//sqared error function
float Picture::error(Picture &image){
  float error=0;
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      error+=pow((get_intensity(j,i)-image.get_intensity(j,i)),2);
    }
  }
  return error;
}

//covariance like error function(the 2nd)
float Picture::error_covariance_like(Picture& image){
  float error=0;
  float divided,divisor_f,divisor_g=0;
  float mean_clean=0;
  float mean_mooved=0;
  int k=0;
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      mean_clean+=get_intensity(j,i);
      mean_mooved+=image.get_intensity(j,i);
      k++;
    }
  }
  mean_clean=mean_clean/k;
  mean_mooved=mean_mooved/k;

  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      divided+=(get_intensity(j,i)-mean_clean)*(image.get_intensity(j,i)-mean_mooved);
      divisor_f+=pow(get_intensity(j,i)-mean_clean,2);
      divisor_g+=pow(image.get_intensity(j,i)-mean_mooved,2);
    }
  }
  error=divided/(sqrt(divisor_f)*sqrt(divisor_g));
  return error;
}


//write in a file.txt the error at each discrete translation.
float Picture::print_loss_function_x_translation(Picture &translated){
  vector <float> errors;
  int x_len= int(x_length);
  int j=0;
  for(int x=(-x_len+1);x<x_len;x++){
    errors.push_back(translation_x(x).error(translated));
  }

  //Now we will write coordinates in a txt in order to plot with python
  ofstream topython("topython_1D.txt", ios::out );

  if(topython){
    int i=0;
    for(int x=(-x_len+1);x<x_len;x++){
      topython<<x<<endl<<errors[i]<<endl;
      i++;
    }
    topython.close();
  }
  else{
    cerr << "Error. Cannot open topython_1D.txt !" << endl;
  }
  //index of the min error.
  return min_element(errors.begin(),errors.end()) - errors.begin() -x_len+1;
}



Point Picture::print_loss_function_xy_translation(Picture &translated){
  vector <float> errors;
  vector<Point> index;
  int j=0;

  for(int x=-(int)x_length+1;x<(int)x_length-1;x++){
    for(int y=-(int)y_length+1;y<(int)y_length-1;y++){

      errors.push_back(translation_opti_int_xy(x,y).error(translated));
      //errors.push_back(floating_translation(x,y).error(translated));

      index.push_back(Point(x,y));
    }
  }


  //Now we will write coordinates in a txt in order to  plot with python
  ofstream topython("topython_2D.txt", ios::out );

  if(topython)
  {
    int i=0;
    for(int x=-(int)x_length+1;x<(int)x_length-1;x++){
      for(int y=-(int)y_length+1;y<(int)y_length-1;y++){
        topython<<x<<endl<<y<<endl<<errors[i]<<endl;
        i++;
      }
    }
    topython.close();
  }
  else{
    cerr << "Error. Cannot open topython_2D.txt !" << endl;
  }
  return index[min_element(errors.begin(),errors.end()) - errors.begin()];
}

//we will translate only arround the difference between 2 center_x coord
//It is an optimization of the previous method
//We should make a confidence interval to make it works with big translations
//which could change the barycenter.
float Picture::loss_function_xt_by_barycenter(Picture &translated){
  vector <float> errors;
  vector<float> index;

  Point pic_center =center_of_pressure();
  Point translated_center=translated.center_of_pressure();

  int interval_center=translated_center.x-pic_center.x;

  int x_len= int(x_length);
  int j=0;

  for(int x=  max( -x_len+1 , interval_center-(x_len/10) )  ; x<min(  x_len,  interval_center+(x_len/10));  x++){
    errors.push_back(translation_x(x).error(translated));
    index.push_back(x);
  }
  //Now we will write coordinates in a txt in order t  plot with python
  ofstream topython("topython_1D_opti.txt", ios::out );

  if(topython){

    int i=0;
    for(int x=  max( -x_len+1 , interval_center-(x_len/10) )  ; x<min(  x_len,  interval_center+(x_len/10));  x++){
      topython<<x<<endl<<errors[i]<<endl;
      i++;
    }

    topython.close();
  }
  else{
    cerr << "Error. Cannot open topython_1D_opti.txt !" << endl;
  }
  return index[min_element(errors.begin(),errors.end()) - errors.begin()];

}

//we will translate only arround the difference between 2 center_x coord and y_center coord
//It is an optimization of the previous method
Point Picture::loss_function_xyt_by_barycenter(Picture &translated){
  vector <float> errors;
  vector<Point> index;

  int x_len= int(x_length);
  int y_len= int(y_length);
  Point pic_center =center_of_pressure();
  Point translated_center=translated.center_of_pressure();
  int interval_center_x=translated_center.x-pic_center.x;
  int interval_center_y=translated_center.y-pic_center.y;

  for(int x=  max( -x_len+1 , interval_center_x-(x_len/10) )  ; x<min(  x_len,  interval_center_x+(x_len/10));  x++){
    for(int y=  max( -y_len+1 , interval_center_y-(y_len/10) )  ; y<min(  y_len,  interval_center_y+(y_len/10));  y++){
      errors.push_back(translation_y(y).translation_x(x).error(translated));
      index.push_back(Point(x,y));

    }
  }
  //Now we will write coordinates in a txt in order to  plot with python
  ofstream topython("topython_2D_opti.txt", ios::out );

  if(topython)
  {
    int i=0;
    for(int x=  max( -x_len+1 , interval_center_x-(x_len/10) )  ; x<min(  x_len,  interval_center_x+(x_len/10));  x++){
      for(int y=  max( -y_len+1 , interval_center_y-(y_len/10) )  ; y<min(  y_len,  interval_center_y+(y_len/10));  y++){
        topython<<x<<endl<<y<<endl<<errors[i]<<endl;
        i++;
      }
    }
    topython.close();
  }
  else{
    cerr << "Error. Cannot open topython_2D_opti.txt !" << endl;
  }
  return index[min_element(errors.begin(),errors.end()) - errors.begin()];

}

Point Picture::loss_function_xyt_by_barycenter_covariance_error(Picture translated){

  vector <float> errors;
  vector<Point> index;
  int x_len= int(x_length);
  int y_len= int(y_length);
  Point pic_center =center_of_pressure();
  Point translated_center=translated.center_of_pressure();
  int interval_center_x=translated_center.x-pic_center.x;
  int interval_center_y=translated_center.y-pic_center.y;

  for(int x=  max( -x_len+1 , interval_center_x-(x_len/10) )  ; x<min(  x_len,  interval_center_x+(x_len/10));  x++){
    for(int y=  max( -y_len+1 , interval_center_y-(y_len/10) )  ; y<min(  y_len,  interval_center_y+(y_len/10));  y++){
      errors.push_back(translation_y(y).translation_x(x).error_covariance_like(translated));
      index.push_back(Point(x,y));
    }
  }
  //Now we will write coordinates in a txt in order to  plot with python
  ofstream topython("topython_2D_opti_covariance.txt", ios::out );

  if(topython)
  {
    int i=0;
    for(int x=  max( -x_len+1 , interval_center_x-(x_len/10) )  ; x<min(  x_len,  interval_center_x+(x_len/10));  x++){
      for(int y=  max( -y_len+1 , interval_center_y-(y_len/10) )  ; y<min(  y_len,  interval_center_y+(y_len/10));  y++){
        topython<<x<<endl<<y<<endl<<errors[i]<<endl;
        i++;
      }
    }
    topython.close();
  }
  else{
    cerr << "Error. Cannot open topython_2D_opti_covariance.txt !" << endl;
  }
  return index[max_element(errors.begin(),errors.end()) - errors.begin()];

}

//this allow us to find with more accuracy the translation coefficient using
//the one obtained with print_xxxxx functions.
float Picture::find_opti_px(float aproxim,Picture &translated){
  vector <float> errors;
  vector <float> nombres_associes;
  int j=0;
  for(int k=1;k<20;k++){
    for(int x=-10;x<10;x++){
      errors.push_back(floating_translation(aproxim+(float)x*pow(10,-k),0).error(translated));
      nombres_associes.push_back(aproxim+(float)x*pow(10,-k));
    }
    int minElementIndex = min_element(errors.begin(),errors.end()) - errors.begin();
    aproxim=nombres_associes[minElementIndex];
    errors.clear();
    nombres_associes.clear();
  }
  return aproxim;
}

//This allow us to find with more accuracy the translation coefficients px and
//py using the ones obtained with print_xxxxx functions.
Point2f Picture::find_opti_px_py(Point_<float> aproxim,Picture &translated){
  vector <float> errors;
  vector <Point2f> index;
  int j=0;
  for(int k=1;k<8;k++){
    for(int x=-10;x<10;x++){
      for(int y=-10;y<10;y++){
        j++;
        errors.push_back(floating_translation(aproxim.x+(float)x*pow(10,-k),aproxim.y+(float)y*pow(10,-k)).error(translated));
        index.push_back(Point2f(aproxim.x+(float)x*pow(10,-k),aproxim.y+(float)y*pow(10,-k)));

      }
    }
    int minElementIndex = min_element(errors.begin(),errors.end()) - errors.begin();
    aproxim=index[minElementIndex];

    errors.clear();
    index.clear();
    j=0;
  }
  return aproxim;
}
