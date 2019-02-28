#include "Picture.h"
#include "Useful_functions.h"


using namespace cv;
using namespace std;

RNG rng(0);

Picture Picture::symmetry_wrt_y()const{
  Picture symmetry;
  symmetry=clone();
  for(int j=0;j<y_length;j++){
    for(int i=0;i<x_length;i++){
      symmetry.set_intensity(j,x_length-1-i,get_intensity(j,i));
    }
  }

  return symmetry;
}


Picture Picture::symmetry_wrt_x()const{
  Picture symmetry;
  symmetry=clone();
  for(int j=0;j<y_length;j++){
    for(int i=0;i<x_length;i++){
      symmetry.set_intensity(y_length-1-j,i,get_intensity(j,i));
    }
  }

  return symmetry;
}
Picture Picture::diagonal_symmetry_top_to_bottom()const{
  Picture sym(picture.t());
  return sym;
}

Picture Picture::diagonal_symmetry_bottom_to_top()const{
  Picture sym(picture.t());
  return sym.symmetry_wrt_y().symmetry_wrt_x();
}
vector<Point> Picture::ellipse_nbh(Point p, unsigned int a, unsigned int b)const{
  int c1 = p.x;
  int c2 = p.y;
  vector<Point> nbh;
  for(unsigned int i=c1-a; i<=c1+a; i++){
    for(unsigned int j=c2-b; j<=c2+b; j++){
      if((double)((i-c1)*(i-c1)/pow(a,2) + (j-c2)*(j-c2)/pow(b,2)) <= 1){
        nbh.push_back(Point(i,j));
      }
    }
  }
  return nbh;
}

void Picture::show_nbh(vector<Point> nbh)const{
  Picture pic_w_nbh;
  pic_w_nbh = this->clone();
  for(Point &p : nbh){
    pic_w_nbh.set_intensity(p.y, p.x, 1);
  }
  pic_w_nbh.print_picture();
}

Picture Picture::transform_isotropic(Point center)const{
  Picture pressure_pic = this->clone();
  for(int x=0;x<x_length;x++){
    for(int y=0;y<y_length;y++){
      Point pixel(x,y);
      float r = norm(center-pixel);

      //float c = 1/(1+0.01*pow(r,1));
      //float c = 1/(1+0.0001  *pow(r,3.));
      //float c = exp(-0.05*pow(r,1));
      float c = exp(-0.0000001*pow(r,4));

      float f = pressure_pic.get_intensity(pixel.y, pixel.x);
      float transform=c*(f-1)+1;
      pressure_pic.set_intensity(pixel.y, pixel.x, transform);
    }
  }
  return pressure_pic;
}

Picture Picture::transform_anisotropic(Point center, unsigned int a, unsigned int b)const{
  Picture pressure_pic = this->clone();
  for(int x=0;x<x_length;x++){
    for(int y=0;y<y_length;y++){
      Point pixel(x,y);
      float r = sqrt((double)((x-center.x)*(x-center.x)/pow(a,2) + (y-center.y)*(y-center.y)/pow(b,2)));
      //float c = 1/(1+1*pow(r,1));
      //float c = 1/(1+1*pow(r,10));
      //float c = exp(-1.5*pow(r,1));
      float c = exp(-1*pow(r,10));
      float f = pressure_pic.get_intensity(pixel.y, pixel.x);
      float transform=c*(f-1)+1;
      pressure_pic.set_intensity(pixel.y, pixel.x, transform);
    }
  }
  return pressure_pic;
}

Picture Picture::extract_ellipse_pic(Point& center, unsigned int a,unsigned int b)const{
  int c1=center.x;
  int c2=center.y;

  Picture extraction=clone();
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      if((double)((i-c1)*(i-c1)/pow(a,2) + (j-c2)*(j-c2)/pow(b,2)) >= 1){
        extraction.set_intensity(j,i,1);
      }
    }
  }
  return extraction;
}

vector<Point> Picture::weak_pressure_border(Point center, unsigned int a, unsigned int b )const{


  Picture pressure_pic(Mat(y_length,x_length,CV_8UC1,Scalar(255)));
  Picture temp_pic(Mat(y_length,x_length,CV_8UC1,Scalar(255)));
  vector<Point> border;
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      double cst = sqrt((double)((i-center.x)*(i-center.x)/pow(a,2) + (j-center.y)*(j-center.y)/pow(b,2)));
      if( cst <= 1 ){
        pressure_pic.set_intensity(j,i,get_intensity(j,i));

      }
      if(cst<=1.001 && cst >=0.999 /*cst==1*/  ){
        if(get_intensity(j,i)<=1){

          border.push_back(Point(i,j));
          }
        }

    }
  }
  vector<Point2f> border_polar_wrt_center;
  vector<Point> new_border;
  vector<Point> sorted_polar;
  vector<float> vector_of_teta;
  vector<Point> random_border;
  for (int c=0;c<border.size();c++){

    float r=norm(border[c]-center);
    float teta=atan2(float(border[c].y-center.y),float(border[c].x-center.x));
    double variance = 3;
    float random_gauss = abs((rng.gaussian(variance))+r/4);
    float new_r=r+random_gauss;
    border_polar_wrt_center.push_back(Point2f(new_r,teta));
    new_border.push_back(Point(int(new_r*cos(teta))+center.x,int(new_r*sin(teta))+center.y));



  }

  sort(border_polar_wrt_center.begin(),border_polar_wrt_center.end(),compare_polar_cord);
  for(int c=0;c<border_polar_wrt_center.size();c++){

    new_border[c].x=border_polar_wrt_center[c].x*cos(border_polar_wrt_center[c].y)+center.x;
    new_border[c].y=border_polar_wrt_center[c].x*sin(border_polar_wrt_center[c].y)+center.y;
    if(new_border[c].x<x_length && new_border[c].x>=0 && new_border[c].y >=0 && new_border[c].y < y_length ){
      pressure_pic.set_intensity(new_border[c].y,new_border[c].x,get_intensity(new_border[c].y,new_border[c].x));
    }
  }

  for(int i=0;i<new_border.size();i++){
    vector<Point> line;
    if(i==new_border.size()-1){
      line =segment(Point(new_border[i]),Point(new_border[0]));
    }
    else{
      line =segment(Point(new_border[i]),Point(new_border[i+1]));
    }
    for(int c=0; c<line.size();c++){
      random_border.push_back(Point(line[c].x,line[c].y));
      if(line[c].x<x_length && line[c].x>=0 && line[c].y >=0 && line[c].y < y_length ){
        pressure_pic.set_intensity(line[c].y,line[c].x,get_intensity(line[c].y,line[c].x));
      }
    }
  }


  return(random_border);
}

vector<Point> Picture::weak_pressure_area(vector<Point> border, Point center, unsigned int a, unsigned int b )const{
  Picture temp_pic(Mat(y_length,x_length,CV_8UC1,Scalar(255)));
  vector<Point> area;
  sort(border.begin(),border.end(),compare_y_cord);
  for(int i=0;i<border.size();i++){
    if(border[i].x<x_length && border[i].x>=0 && border[i].y >=0 && border[i].y < y_length ){
    temp_pic.set_intensity(border[i].y,border[i].x,0);
  }
    area.push_back(border[i]);
  }


  for(int i=0;i<border.size()-1;i++){
      if(temp_pic.get_intensity(border[i].y-1,border[i].x+1) == 0 && border[i].y == border[i+1].y ){
      vector<Point> line;
      line = segment(border[i],border[i+1]);
      for(int j=0;j<line.size();j++){
        double is_elipse=(double)((line[j].x-center.x)*(line[j].x-center.x)/pow(a,2) + (line[j].y-center.y)*(line[j].y-center.y)/pow(b,2));
        if(is_elipse >= 1){
          area.push_back(Point(line[j].x,line[j].y));
          temp_pic.set_intensity(line[j].y,line[j].x,0);
        }
      }
    }
  }
  return area;
}


Picture Picture::attenuation_weak_area(Point center, unsigned int a, unsigned int b)const{
  vector<Point> border_area = weak_pressure_border(center,a,b);
  vector<Point> area = weak_pressure_area(border_area,center,a,b);
  Picture pressure_pic(Mat(y_length,x_length,CV_8UC1,Scalar(255)));
  vector<Point> Elipse = ellipse_nbh(center,a,b);
  vector<Point> Ellipse_border;
  float min_distance = min(a,b);
  float max_distance = max(a,b);
  // Get the border of the Ellipse
  for (int i=0;i<Elipse.size();i++){
    /*if(get_intensity(Elipse[i].y,Elipse[i].x)>0.5){
      pressure_pic.set_intensity(Elipse[i].y,Elipse[i].x,1);
    }
    else{
      */pressure_pic.set_intensity(Elipse[i].y,Elipse[i].x,get_intensity(Elipse[i].y,Elipse[i].x));
    //}*/
    double cst = (double)((Elipse[i].x-center.x)*(Elipse[i].x-center.x)/pow(a,2) + (Elipse[i].y-center.y)*(Elipse[i].y-center.y)/pow(b,2));

    if(cst<1.1 && cst>0.9){
      Ellipse_border.push_back(Elipse[i]);
    }
  }



  for (int i=0; i<area.size();i++){
    // We compute the distance between points in the random area and the Ellipse.
    float min_distance_ellipse=norm(area[i]-Ellipse_border[0]);
    pressure_pic.set_intensity(area[i].y,area[i].x,1);
    for(int j=1; j<Ellipse_border.size();j++){
      float Norm_el=norm(area[i]-Ellipse_border[j]);
      if(Norm_el<min_distance_ellipse){
        min_distance_ellipse=Norm_el;
      }
    }


    // We compute the distance between points in the random area and the "border of the random area"
    float min_distance_border=norm(area[i]-border_area[0]);
    for(int j=1; j<border_area.size();j++){
      float Norm_bord=norm(area[i]-border_area[j]);
      pressure_pic.set_intensity(border_area[j].y,border_area[j].x,1);
      if(Norm_bord<min_distance_border){
        min_distance_border=Norm_bord;
      }
    }
    float intensity = get_intensity(area[i].y,area[i].x);
    float normalisation_factor=1/(pow(pow(min_distance_border,1.8)+pow(min_distance_ellipse,1),1  ));
    float strange_distance=pow(min_distance_ellipse,1)*normalisation_factor;
    //strange_distance=strange_distance;
    //float transformation=1-(1-intensity)*pow(strange_distance,4);

    //float transformation_distance=(exp(-pow(1-strange_distance,1))-1)/(exp(-1)-1);
    float transformation_distance=pow((1-strange_distance),1);
    float transformation_intensity=intensity;
    //float transformation=strange_distance*(intensity-1)+1;
    //float transformation=strange_distance;
    float transformation=transformation_distance*(transformation_intensity-1)+1;
    //pressure_pic.set_intensity(area[i].y,area[i].x,get_intensity(area[i].y,area[i].x));
    pressure_pic.set_intensity(area[i].y,area[i].x,transformation);
    /*if(pressure_pic.get_intensity(area[i].y,area[i].x)>0.5 && pressure_pic.get_intensity(area[i].y+1,area[i].x)>0.2 && pressure_pic.get_intensity(area[i].y-1,area[i].x)>0.2 &&  pressure_pic.get_intensity(area[i].y,area[i].x+1)>0.2 &&  pressure_pic.get_intensity(area[i].y,area[i].x-1)>0.2){
      pressure_pic.set_intensity(area[i].y,area[i].x,1);
    }*/
  }
  return(pressure_pic);
}
