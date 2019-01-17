#include <stdio.h>
#include <iostream>
#include "Picture.h"

using namespace cv;
using namespace std;


int main(int argc, char** argv )
{
  float temps;
  clock_t t0,t1, t2,t_end;
  t0 = clock();

  Picture clean("/home/tristan/Bureau/Projet/Digital-Analysis-of-Fingerprints-MSIAM/infra/data/clean_finger.png");
  Picture tx_finger("/home/tristan/Bureau/Projet/Digital-Analysis-of-Fingerprints-MSIAM/infra/data/tx_finger.png");
  Picture txy_finger("/home/tristan/Bureau/Projet/Digital-Analysis-of-Fingerprints-MSIAM/infra/data/txy_finger.png");


  //Draw the loss function with respect to px between clean and tx_finger.
  t1=clock();
  float px_brute_force=clean.print_loss_function_x_translation(tx_finger);
  t2=clock();
  cout<<"With brute_force method we obtain px = "<<px_brute_force<<endl;
  cout<<"It took "<<(float)(t2-t1)/CLOCKS_PER_SEC<<" seconds to run."<<endl;

  cout<<"--------------------------------------------------------"<<endl;
  //Draw the loss function with respect to px & py between clean and txy_finger.
  //Brute force method, it take between 9 and 15 min to run depending on
  //the translation method you use

  t1=clock();
  Point Pxy_brute_force=clean.print_loss_function_xy_translation(txy_finger);
  t2=clock();
  cout<<"With brute_force method we obtain px = "<<Pxy_brute_force.x;
  cout<<" py = "<<Pxy_brute_force.y<<endl;
  cout<<"It took "<<(float)(t2-t1)/CLOCKS_PER_SEC<<" seconds to run."<<endl;

  cout<<"--------------------------------------------------------"<<endl;
  //Draw the loss function with respect to px between clean and tx_finger.
  //Smarter method.
  t1=clock();
  float px_smart = clean.loss_function_xt_by_barycenter(tx_finger);
  t2=clock();
  cout<<"With smarter method we obtain px = "<<px_smart<<endl;
  cout<<"It took "<<(float)(t2-t1)/CLOCKS_PER_SEC<<" seconds to run."<<endl;

  cout<<"--------------------------------------------------------"<<endl;
  //Draw the loss function with respect to px and py between clean and txy_finger.
  //Smarter method, it doesn't take more than 1 min to run with the non
  //optimized translation method.
  t1=clock();
  Point Pxy_smart=clean.loss_function_xyt_by_barycenter(txy_finger);
  t2=clock();
  cout<<"With smarter method we obtain px = "<<Pxy_smart.x;
  cout<<" py = "<<Pxy_smart.y<<endl;
  cout<<"It took "<<(float)(t2-t1)/CLOCKS_PER_SEC<<" seconds to run."<<endl;


  cout<<"--------------------------------------------------------"<<endl;
  //Draw the 2nd loss function with respect to px and py between clean and
  //txy_finger with the fastest method we have.
  t1=clock();
  Point Pxy_smart_2nd_error=clean.loss_function_xyt_by_barycenter_covariance_error(txy_finger);
  t2=clock();
  cout<<"With Smarter method  and 2nd error functions we obtain px = "<<Pxy_smart_2nd_error.x;
  cout<<" py = "<<Pxy_smart_2nd_error.y<<endl;
  cout<<"It took "<<(float)(t2-t1)/CLOCKS_PER_SEC<<" seconds to run."<<endl;


  cout<<"--------------------------------------------------------"<<endl;
  //The absolute error image ( saved as "absolute_error_image_int_px_py.png")
  Picture abs_error=txy_finger-clean.floating_translation(Pxy_smart_2nd_error.x,Pxy_smart_2nd_error.y);
  abs_error.print_picture();
  abs_error.SAVE_PIC("absolute_error_image_int_px_py.png");


  cout<<"--------------------------------------------------------"<<endl;
  //To find exact translation coefficient for x translation only :
  float exact_px_coeff=clean.find_opti_px(px_smart,tx_finger);
  cout<<"exact px : "<<exact_px_coeff<<endl;

  cout<<"--------------------------------------------------------"<<endl;
  //To find exact translation coefficient for x&y translation :
  Point2f exact_px_py_coeff=clean.find_opti_px_py(Pxy_smart,txy_finger);
  cout<<"exact px and py coeffs : "<<exact_px_py_coeff<<endl;


  cout<<"--------------------------------------------------------"<<endl;
  //The absolute error image with optimized px and py
  // saved as "absolute_error_image_float_px_py.png"
  Picture abs_error_opti=txy_finger-clean.floating_translation(exact_px_py_coeff.x,exact_px_py_coeff.y);
  abs_error_opti.print_picture();
  abs_error_opti.SAVE_PIC("absolute_error_image_float_px_py.png");


  cout<<"--------------------------------------------------------"<<endl;


  t_end = clock();
  temps = (float)(t_end-t0)/CLOCKS_PER_SEC;
  cout<<"run time : "<<temps<<endl;
  return 0;
}
