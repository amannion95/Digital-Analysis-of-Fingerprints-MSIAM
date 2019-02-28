#include "Picture.h"
#include "Useful_functions.h"



using namespace cv;
using namespace std;
using namespace Eigen;

Picture Picture::discrete_convolution(Matrix<float,Dynamic,Dynamic> mask)const{
  Matrix<float,Dynamic,Dynamic> image;
  cv2eigen(picture,image);
  int MaskCenterX=((mask.cols())/2);
  int MaskCenterY=((mask.rows())/2);
  Matrix<float,Dynamic,Dynamic> convolution(y_length,x_length);
  for(int i=0; i<y_length;i++){
    for(int j=0; j<x_length;j++){
      convolution(i,j)=0;
      for(int m=0;m<mask.rows();m++){
        int mm=mask.rows()-1-m;
        for(int n=0;n<mask.cols();n++){
          int nn=mask.cols()-1-n;
          int ii = i+(m-MaskCenterY);
          int jj = j+(n-MaskCenterX);
          if(ii<0){
            ii++;
          }
          if(jj<0){
            jj++;
          }
          if(ii>=image.rows()){
            ii=ii-1;
          }
          if(jj>=image.cols()){
            jj=jj-1;
          }
          if(ii>=0 && ii<y_length && jj>=0 && jj<x_length){
            convolution(i,j)= convolution(i,j)+float(image(ii,jj)*mask(mm,n));
          }

        }
      }
    }
  }
  Mat convolution_dst, convolution_dst_8UCV1;
  eigen2cv(convolution,convolution_dst);
  convolution_dst.convertTo(convolution_dst_8UCV1,CV_8UC1);
  return(Picture(convolution_dst_8UCV1));
}
/*
Picture Picture::ConvolutionDFT(Matrix<float,Dynamic,Dynamic> mask){


// Picture est notre classe pour manipuler des images, picture ( de Type Mat ) en est un attribut qui stock la matrice de l'image.
// Le mask est donné dans une matrice dynamique de eigen.
  Mat picture_temp,mask_cv,convolution,picture_temp_dst,mask_cv_temp;


  eigen2cv(mask,mask_cv); // On converti le mask en objet MaT
  picture.convertTo(picture_temp,CV_32FC1); // dft ne marche pas avec n'importe quel type mat (pas avec CV_8UC1). on convertis donc tout en CV_32FC1
  mask_cv.convertTo(mask_cv_temp,CV_32FC1);     // idem
  dft(picture_temp,picture_temp,0,0); // on fait notre dft de picture_temp qu'on stock dans picture_temp
  dft(mask_cv_temp,mask_cv_temp,0,0);  // idem
  copyMakeBorder(mask_cv_temp,mask_cv_temp,0,x_length-mask.rows(),0,y_length-mask.cols(),BORDER_CONSTANT,0.); // on ajoute des 0 au mask pour rendre la multiplication matriciel possible
  convolution=picture_temp*mask_cv_temp;  // on fait le produit matriciel
  idft(convolution,convolution); // On applique la dft inverse
  convolution.convertTo(convolution,CV_8UC1); // on reconvertis la convolution en CV_8UC1


  return(Picture(convolution));

}*/

Picture Picture::ConvolutionDFT(Matrix<float,Dynamic,Dynamic> mask)const{
  Mat src, mask_src, convolution;
  picture.convertTo(src,CV_64FC1);
  eigen2cv(mask,mask_src);
  mask_src.convertTo(mask_src,CV_64FC1);
  convolveDFT(src,mask_src,convolution);
  convolution.convertTo(convolution,CV_8UC1);
  cout << convolution.rows << " and " << convolution.cols << "mask size : " << mask.rows() << '\n';
  int n = mask.rows()/2;
  Rect myRegion(/*mask.cols()/2*/n,/*mask.rows()/2*/n,x_length,y_length);
  Mat convolution_reshaped;
  convolution_reshaped = convolution(myRegion);

  return(Picture(convolution_reshaped));
}

Picture Picture::discrete_convolution_evolutive_kernel_no_identity()const{
  Matrix<float,Dynamic,Dynamic> image;
  cv2eigen(picture,image);


  Matrix<float,Dynamic,Dynamic> convolution(y_length,x_length);
  float size = 7;
  Point center = center_of_pressure();
  float max_distance_x=max(center.x,int(x_length-center.x));
  float max_distance_y=max(center.y,int(y_length-center.y));
  //float max_distance=max(max_distance_x,max_distance_y);
  float max_distance=1 ;
  float a=125;
  float b=150;
  for(int i=0; i<y_length;i++){
    for(int j=0; j<x_length;j++){
      convolution(i,j)=0;
      float distance = pow(j-center.x,2)/pow(a,2) + pow(i-center.y,2)/pow(b,2);
      //float distance=norm(center-Point(j,i));
      float normalized_distance=distance/max_distance;
      Matrix<float,Dynamic,Dynamic> mask(Evolutive_kernel_no_identity(normalized_distance,size));
      int MaskCenterX=((mask.cols())/2);
      int MaskCenterY=((mask.rows())/2);
      if(distance<1){
      for(int m=0;m<mask.rows();m++){
        int mm=mask.rows()-1-m;
        for(int n=0;n<mask.cols();n++){
          int nn=mask.cols()-1-n;
          int ii = i+(m-MaskCenterY);
          int jj = j+(n-MaskCenterX);
          if(ii<0){
            ii++;
          }
          if(jj<0){
            jj++;
          }
          if(ii>=image.rows()){
            ii=ii-1;
          }
          if(jj>=image.cols()){
            jj=jj-1;
          }
          if(ii>=0 && ii<y_length && jj>=0 && jj<x_length){
            convolution(i,j)= convolution(i,j)+float(image(ii,jj)*mask(mm,n));
          }

        }
      }
    }
    else(convolution(i,j)=255);
    }
  }
  Mat convolution_dst, convolution_dst_8UCV1;
  eigen2cv(convolution,convolution_dst);
  convolution_dst.convertTo(convolution_dst_8UCV1,CV_8UC1);
  return(Picture(convolution_dst_8UCV1));
}

Picture Picture::discrete_convolution_evolutive_kernel()const{
  Matrix<float,Dynamic,Dynamic> image;
  cv2eigen(picture,image);


  Matrix<float,Dynamic,Dynamic> convolution(y_length,x_length);
  float size = 7;
  Point center = center_of_pressure();
  float max_distance_x=max(center.x,int(x_length-center.x));
  float max_distance_y=max(center.y,int(y_length-center.y));
  //float max_distance=max(max_distance_x,max_distance_y);
  float max_distance=1 ;
  float a=125;
  float b=150;
  for(int i=0; i<y_length;i++){
    for(int j=0; j<x_length;j++){
      convolution(i,j)=0;
      float distance = pow(j-center.x,2)/pow(a,2) + pow(i-center.y,2)/pow(b,2);
      //float distance=norm(center-Point(j,i));
      float normalized_distance=distance/max_distance;
      Matrix<float,Dynamic,Dynamic> mask(Evolutive_kernel(normalized_distance,size));
      int MaskCenterX=((mask.cols())/2);
      int MaskCenterY=((mask.rows())/2);
      if(distance<1){
      for(int m=0;m<mask.rows();m++){
        int mm=mask.rows()-1-m;
        for(int n=0;n<mask.cols();n++){
          int nn=mask.cols()-1-n;
          int ii = i+(m-MaskCenterY);
          int jj = j+(n-MaskCenterX);
          if(ii<0){
            ii++;
          }
          if(jj<0){
            jj++;
          }
          if(ii>=image.rows()){
            ii=ii-1;
          }
          if(jj>=image.cols()){
            jj=jj-1;
          }
          if(ii>=0 && ii<y_length && jj>=0 && jj<x_length){
            convolution(i,j)= convolution(i,j)+float(image(ii,jj)*mask(mm,n));
          }

        }
      }
    }
    else(convolution(i,j)=255);
    }
  }
  Mat convolution_dst, convolution_dst_8UCV1;
  eigen2cv(convolution,convolution_dst);
  convolution_dst.convertTo(convolution_dst_8UCV1,CV_8UC1);
  return(Picture(convolution_dst_8UCV1));
}
