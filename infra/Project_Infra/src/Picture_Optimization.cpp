#include "Picture.h"
#include "Usefull_functions.h"

#include <time.h>
#include <fstream>

using namespace cv;
using namespace std;


Picture Picture::translation_x(int coeff)const{
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

Picture Picture::translation_y(int coeff)const{
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


Picture Picture::translation_opti_int_xy(int x,int y)const{
  int offset_x=max(0,-x);
  int offset_y=max(0,-y);

  Rect roi;

  //We choose as region of interest the part of the picture which will still
  //be in the picture after the translation
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
Picture Picture::floating_translation(float tx, float ty)const{
  Mat image2;
  Mat M=(Mat_<float>(2,3)<<1,0,tx,0,1,ty);
  warpAffine( picture, image2, M, picture.size());
  Picture translated(image2);

  //warpAffine is pading with black pixels instead of white pixels.
  //These loops are correcting it
  if (tx>0){
    for(int x=0;x<tx;x++){
      for(int y=0;y<y_length;y++){
        translated.set_intensity(y,x,1);
      }
    }
  }
  else{
    for(int x=x_length-1;x>x_length-1+tx;x--){
      for(int y=0;y<y_length;y++){
        translated.set_intensity(y,x,1);
      }
    }
  }

  if (ty>0){
    for(int y=0;y<ty;y++){
      for(int x=0;x<x_length;x++){
        translated.set_intensity(y,x,1);
      }
    }
  }
  else{
    for(int y=y_length-1;y>y_length-1+ty;y--){
      for(int x=0;x<x_length;x++){
        translated.set_intensity(y,x,1);
      }
    }
  }
  return translated;
}


bool Picture::is_same(Picture& image)const{
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      if (get_intensity(j,i)!=image.get_intensity(j,i)){
        return false;
      }
    }
  }
  return true;
}

Picture Picture::operator-(Picture& to_substract)const{
  float ab;

  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      ab=abs(get_intensity(j,i)-to_substract.get_intensity(j,i));
      to_substract.set_intensity(j,i,ab);
    }
  }
  return to_substract;
}

//sum of sqared error function
float Picture::error(Picture &image)const{
  float error=0;
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      error+=pow((get_intensity(j,i)-image.get_intensity(j,i)),2);
    }
  }
  return error;
}

//covariance like error function(the 2nd)
float Picture::error_covariance_like(Picture& image)const{
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
      divided+=(get_intensity(j,i)-mean_clean)*(image.get_intensity(j,i)\
                                                -mean_mooved);
      divisor_f+=pow(get_intensity(j,i)-mean_clean,2);
      divisor_g+=pow(image.get_intensity(j,i)-mean_mooved,2);
    }
  }
  error=divided/(sqrt(divisor_f)*sqrt(divisor_g));
  return error;
}


//write in a file.txt the error at each discrete translation.
float Picture::print_loss_function_x_translation(Picture &translated)const{
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



Point Picture::print_loss_function_xy_translation(Picture &translated)const{
  vector <float> errors;
  vector<Point> index;
  int j=0;

  for(int x=-(int)x_length+1;x<(int)x_length-1;x++){
    for(int y=-(int)y_length+1;y<(int)y_length-1;y++){
      errors.push_back(translation_opti_int_xy(x,y).error(translated));
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

float Picture::loss_function_xt_by_barycenter(Picture &translated)const{
  vector <float> errors;
  vector<float> index;

  Point pic_center =center_of_pressure();
  Point translated_center=translated.center_of_pressure();

  int interval_center=translated_center.x-pic_center.x;

  int x_len= int(x_length);
  int j=0;

  for(int x=  max( -x_len+1 , interval_center-(x_len/10) );
          x<min(  x_len,  interval_center+(x_len/10));  x++){
    errors.push_back(translation_x(x).error(translated));
    index.push_back(x);
  }
  //Now we will write coordinates in a txt in order t  plot with python
  ofstream topython("topython_1D_opti.txt", ios::out );

  if(topython){

    int i=0;
    for(int x=  max( -x_len+1 , interval_center-(x_len/10) )  ;
        x<min(  x_len,  interval_center+(x_len/10));  x++){
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

Point Picture::loss_function_xyt_by_barycenter(Picture &translated)const{
  vector <float> errors;
  vector<Point> index;

  int x_len= int(x_length);
  int y_len= int(y_length);
  Point pic_center =center_of_pressure();
  Point translated_center=translated.center_of_pressure();
  int interval_center_x=translated_center.x-pic_center.x;
  int interval_center_y=translated_center.y-pic_center.y;

  for(int x=  max( -x_len+1 , interval_center_x-(x_len/10) )  ;
      x<min(  x_len,  interval_center_x+(x_len/10));  x++){
    for(int y=  max( -y_len+1 , interval_center_y-(y_len/10) )  ;
        y<min(  y_len,  interval_center_y+(y_len/10));  y++){

      errors.push_back(translation_y(y).translation_x(x).error(translated));
      index.push_back(Point(x,y));
    }
  }
  //Now we will write coordinates in a txt in order to  plot with python
  ofstream topython("topython_2D_opti.txt", ios::out );

  if(topython)
  {
    int i=0;
    for(int x=  max( -x_len+1 , interval_center_x-(x_len/10) )  ;
        x<min(  x_len,  interval_center_x+(x_len/10));  x++){
      for(int y=  max( -y_len+1 , interval_center_y-(y_len/10) )  ;
          y<min(  y_len,  interval_center_y+(y_len/10));  y++){
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

Point Picture::loss_function_xyt_by_barycenter_covariance_error(Picture& translated)const{

  vector <float> errors;
  vector<Point> index;
  int x_len= int(x_length);
  int y_len= int(y_length);
  Point pic_center =center_of_pressure();
  Point translated_center=translated.center_of_pressure();
  int interval_center_x=translated_center.x-pic_center.x;
  int interval_center_y=translated_center.y-pic_center.y;

  for(int x=  max( -x_len+1 , interval_center_x-(x_len/10) )  ;
      x<min(  x_len,  interval_center_x+(x_len/10));  x++){
    for(int y=  max( -y_len+1 , interval_center_y-(y_len/10) )  ;
      y<min(  y_len,  interval_center_y+(y_len/10));  y++){
      errors.push_back(translation_y(y).translation_x(x).error_covariance_like(translated));
      index.push_back(Point(x,y));
    }
  }
  //Now we will write coordinates in a txt in order to  plot with python
  ofstream topython("topython_2D_opti_covariance.txt", ios::out );

  if(topython)
  {
    int i=0;
    for(int x=  max( -x_len+1 , interval_center_x-(x_len/10) )  ;
        x<min(  x_len,  interval_center_x+(x_len/10));  x++){
      for(int y=  max( -y_len+1 , interval_center_y-(y_len/10) )  ;
          y<min(  y_len,  interval_center_y+(y_len/10));  y++){
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

float Picture::find_opti_px(float aproxim,Picture &translated)const{
  vector <float> errors;
  vector <float> nombres_associes;
  int j=0;
  for(int k=1;k<20;k++){
    if (floating_translation(aproxim,0).is_same(translated)){
      return aproxim;
    }
    for(int x=-10;x<10;x++){
      errors.push_back(floating_translation(aproxim+(float)x*pow(10,-k),0)\
                                                        .error(translated));
      nombres_associes.push_back(aproxim+(float)x*pow(10,-k));
    }
    int minElementIndex = min_element(errors.begin(),errors.end()) - errors.begin();
    aproxim=nombres_associes[minElementIndex];
    errors.clear();
    nombres_associes.clear();
  }
  return aproxim;
}

Point2f Picture::find_opti_px_py(Point_<float> aproxim,Picture &translated)const{
  vector <float> errors;
  vector <Point2f> index;
  int j=0;
  for(int k=1;k<8;k++){
    if (floating_translation(aproxim.x,aproxim.y).is_same(translated)){
      cout<<"Ok : "<<aproxim<<endl;
      return aproxim;
    }
    for(int x=-10;x<10;x++){
      for(int y=-10;y<10;y++){
        j++;
        errors.push_back(floating_translation(aproxim.x+(float)x*pow(10,-k),\
                          aproxim.y+(float)y*pow(10,-k)).error(translated));
        index.push_back(Point2f(aproxim.x+(float)x*pow(10,-k)\
                        ,aproxim.y+(float)y*pow(10,-k)));
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

Point2f Picture::estimation_tranlsation_by_dft(Picture& translated)const{
    Mat MAT,MAT2;
    picture.convertTo(MAT, CV_64FC1);
    translated.picture.convertTo(MAT2, CV_64FC1);
    Picture a1 (MAT);
    Picture a2 (MAT2);
    return phaseCorrelate(a1.picture,a2.picture);
}

float Picture::estimation_rotation_bruteforce(Picture& rotated)const{
  vector <float> errors;
  Point middle(x_length/2,y_length/2);
  for(int angle=0;angle<360;angle++){
    errors.push_back(bilinear_rotation_polar(middle,angle).error(rotated));
  }
  return min_element(errors.begin(),errors.end()) - errors.begin();
}


float Picture::more_accurate_rotation_parameter(float angle, Picture& rotated)const{
  vector <float> errors;
  vector <float> nombres_associes;

  Point middle(x_length/2,y_length/2);
  for(int k=1;k<6;k++){
    for(int x=-10;x<10;x++){
      errors.push_back(bilinear_rotation_polar(middle,angle+(float)x*pow(10,-k))\
                      .error(rotated));
      nombres_associes.push_back(angle+(float)x*pow(10,-k));
    }

    int minElementIndex = min_element(errors.begin(),errors.end()) - errors.begin();
    angle=nombres_associes[minElementIndex];
    errors.clear();
    nombres_associes.clear();
  }
  return angle;
}

Picture Picture::put_barycenter_at_picture_center()const{
  Point pic_barycenter =center_of_pressure();
  Point middle(x_length/2,y_length/2);
  return floating_translation(middle.x-pic_barycenter.x,middle.y-pic_barycenter.y);
}

float Picture::sum_intensity_picture()const{
  float res;
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      res+=get_intensity(j,i);
    }
  }
  return res;
}

void Picture::execution_evaluation_rtxy(Picture& Rot_txy_picture)const{

  Point milieu(x_length/2,y_length/2);

  Picture clean_mooved_in_center=put_barycenter_at_picture_center();
  Picture Rtxy_mooved_in_center=Rot_txy_picture.put_barycenter_at_picture_center();

  float sum_clean=clean_mooved_in_center.sum_intensity_picture();
  float sum_rotated=Rtxy_mooved_in_center.sum_intensity_picture();

  cout<<sum_clean<<endl;
  cout<<sum_rotated<<endl;
  cout<<"ratio : "<<sum_clean/sum_rotated<<endl;

  Picture center_blured=Rtxy_mooved_in_center.apply_gaussian_blur(25);

  float angle=clean_mooved_in_center.apply_gaussian_blur(25)\
              .estimation_rotation_bruteforce(center_blured);
  Point2f blur_estim_txy=bilinear_rotation_polar(milieu,angle)\
              .estimation_tranlsation_by_dft(Rot_txy_picture);


  float angle_w_o_blur=clean_mooved_in_center.estimation_rotation_bruteforce\
                                              (Rtxy_mooved_in_center);
  Point2f w_o_blur_txy=bilinear_rotation_polar(milieu,angle_w_o_blur)\
                        .estimation_tranlsation_by_dft(Rot_txy_picture);

  if (bilinear_rotation_polar(milieu,angle_w_o_blur).floating_translation\
          (w_o_blur_txy.x,w_o_blur_txy.y).error(Rot_txy_picture)
      <bilinear_rotation_polar(milieu,angle).floating_translation
          (blur_estim_txy.x,blur_estim_txy.y).error(Rot_txy_picture)){

      cout<<"rot estim "<<angle_w_o_blur<<endl;
      cout<<"txy estimated"<<w_o_blur_txy<<endl;
  }
  else{
    cout<<"rot estim "<<angle<<endl;
    cout<<"txy estimated"<<blur_estim_txy<<endl;
  }
}

//----------------------------Gradient descent functions-----------------------

float** Picture::dg_dwx()const{

  float **dg_px = new float*[y_length];
  for ( int j = 0 ; j < y_length ; j++ ) {
    dg_px[j] = new float[x_length];
  }

  for (int i=1;i<x_length-1;i++){
    for(int j=0;j<y_length;j++){
      dg_px[j][i]=(get_intensity(j,i+1)-get_intensity(j,i));
    }
  }
  for (int i=0;i<y_length;i++){
    dg_px[i][x_length-1]=0;
    dg_px[i][0]=0;
  }
  return dg_px;
}

float** Picture::dg_dwy()const{

  float **dg_py = new float*[y_length];
  for ( int j = 0 ; j < y_length ; j++ ) {
    dg_py[j] = new float[x_length];
  }

  for (int i=0;i<x_length;i++){
    for(int j=1;j<y_length-1;j++){

      dg_py[j][i]=(get_intensity(j+1,i) - get_intensity(j,i));
    }
  }
  for (int i=0;i<x_length;i++){
    dg_py[y_length-1][i]=0;
    dg_py[0][i]=0;

  }
  return dg_py;
}

Point2f Picture::Error_partial_deriv_px_py(Point2f p_after,Picture& translated)const{

  Picture translated_p_after=floating_translation(p_after.x,p_after.y);
  Picture translated_p_after_x=floating_translation(p_after.x,0);
  Picture translated_p_after_y=floating_translation(0,p_after.y);

  float **dg_px = new float*[y_length];
  for ( int j = 0 ; j < y_length ; j++ ) {
    dg_px[j] = new float[x_length];
  }
  float **dg_py = new float*[y_length];
  for ( int j = 0 ; j < y_length ; j++ ) {
    dg_py[j] = new float[x_length];
  }

  dg_px=translated_p_after_x.dg_dwx();
  dg_py=translated_p_after_y.dg_dwy();
  float sum_x=0;
  float sum_y=0;
  for (int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      sum_x=sum_x-2*dg_px[j][i]*(translated_p_after_x\
            .get_intensity(j,i)-translated.get_intensity(j,i));
      sum_y=sum_y-2*dg_py[j][i]*(translated_p_after_y\
            .get_intensity(j,i)-translated.get_intensity(j,i));
    }
  }

  for(int j = 0; j < y_length; j++){
    delete[] dg_px[j];
    delete[] dg_py[j];
  }
  delete[] dg_px;
  delete[] dg_py;

  return Point2f(sum_x,sum_y);
}

Point2f Picture::gradient_descent(Point2f p,float epsilon,float alpha,Picture& translated)const{

  Point2f grad = Error_partial_deriv_px_py(p,translated);

  if ((abs(grad.x)<epsilon)&&(abs(grad.y)<epsilon)){
    return p;
  }

  p.x=p.x-alpha*grad.x;
  p.y=p.y-alpha*grad.y;

  cout<<" p : "<<p<<endl;
  cout<<"grad : "<<grad<<endl;
  cout<<"loss function (should decrease at each step) : "\
      <<floating_translation(p.x,p.y).error(translated)<<endl;
  cout<<"----------------------------------"<<endl;

  return gradient_descent(p,epsilon,alpha,translated);
}
