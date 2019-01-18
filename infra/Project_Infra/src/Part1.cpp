#include "Picture.h"
#include "Usefull_functions.h"


using namespace cv;
using namespace std;

RNG rng(0);

Picture Picture::symmetry_wrt_y()const{
  Picture symmetry;
  symmetry=clone();
  std::cout<<symmetry.x_length<<std::endl;
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
vector<Point> Picture::ellipse_nbh(Point p, unsigned int a, unsigned int b){
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

Picture Picture::log_transform_isotropic(Point p, unsigned int a, unsigned int b, double coef){
  Picture pressure_pic = this->clone();
  vector<Point> nbh = ellipse_nbh(p, a, b);
  for(Point &pixel : nbh){
    float c = log_coeff_isotropic(pixel, p, coef);
    float m = pressure_pic.get_intensity(pixel.y, pixel.x);
    pressure_pic.set_intensity(pixel.y, pixel.x, c*m);
  }
  return pressure_pic;
}

Picture Picture::pow_transform_isotropic(Point p, unsigned int a, unsigned int b, double coef){
  Picture pressure_pic = this->clone();
  vector<Point> nbh = ellipse_nbh(p, a, b);
  for(Point &pixel : nbh){
    float r = norm(p-pixel);
    float m = pressure_pic.get_intensity(pixel.y, pixel.x);
    pressure_pic.set_intensity(pixel.y, pixel.x, pow(1/1.01,r)*m);
  }
  return pressure_pic;
}

Picture Picture::extract_ellipse_pic(Point center, unsigned int a,unsigned int b){
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

Picture Picture::accentuation_diff(int winsize ){
  float mean;
  Picture filtered=clone();
  int semi_size=(winsize-1)/2;
  vector<float> conteneur;
  for(int i=0;i<x_length-semi_size;i++){
    for(int j=0;j<y_length-semi_size;j++){
      for(int u=-semi_size;u<=semi_size;u++){
        for(int v=-semi_size;v<semi_size;v++){
          if ((u!=0) || (v!=0 )){
            conteneur.push_back(get_intensity(j,i)-get_intensity(v,u));
          }
        }
      }

      for(std::vector<float>::iterator it = conteneur.begin(); it != conteneur.end(); ++it){
        mean += *it;
      }
      mean=mean/conteneur.size();
      cout<<"taille : "<<conteneur.size()<<endl;
      if (get_intensity(j,i)+mean>1){
        filtered.set_intensity(j,i,1);
      }
      else if ((get_intensity(j,i)+mean)<0){
        filtered.set_intensity(j,i,0);

      }
      else{
        filtered.set_intensity(j,i,get_intensity(j,i)+mean);
      }
      conteneur.clear();
    }
  }
  return filtered;
}

Picture Picture::without_noise(){
  Picture filtered=clone();

  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      filtered.set_intensity(j,i,min((float)(1+pow(get_intensity(j,i),1))*get_intensity(j,i),(float)1));

    }
  }
  return filtered;
}

vector<Point> Picture::weak_pressure_border(Point center, unsigned int a, unsigned int b ){


  Picture pressure_pic(Mat(y_length,x_length,CV_8UC1,Scalar(255)));
  vector<Point> border;
  for(int i=0;i<x_length;i++){
    for(int j=0;j<y_length;j++){
      double cst = (double)((i-center.x)*(i-center.x)/pow(a,2) + (j-center.y)*(j-center.y)/pow(b,2));
      if( cst <= 1 ){
        pressure_pic.set_intensity(j,i,get_intensity(j,i));
      }
      if(cst<=1.001 && cst >=0.999){
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
    double variance = 1;
    float random_gauss = abs(rng.gaussian(variance)+r/10);
    float new_r=r+random_gauss;
    border_polar_wrt_center.push_back(Point2f(new_r,teta));
    new_border.push_back(Point(int(new_r*cos(teta))+center.x,int(new_r*sin(teta))+center.y));


  }
  sort(border_polar_wrt_center.begin(),border_polar_wrt_center.end(),compare_polar_cord);
  for(int c=0;c<border_polar_wrt_center.size();c++){

    new_border[c].x=border_polar_wrt_center[c].x*cos(border_polar_wrt_center[c].y)+center.x;
    new_border[c].y=border_polar_wrt_center[c].x*sin(border_polar_wrt_center[c].y)+center.y;
    pressure_pic.set_intensity(new_border[c].y,new_border[c].x,get_intensity(new_border[c].y,new_border[c].x));

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
      pressure_pic.set_intensity(line[c].y,line[c].x,get_intensity(line[c].y,line[c].x));
    }
  }

  return(random_border);
}

vector<Point> Picture::weak_pressure_area(vector<Point> border, Point center, unsigned int a, unsigned int b ){
  Picture temp_pic(Mat(y_length,x_length,CV_8UC1,Scalar(255)));
  vector<Point> area;
  sort(border.begin(),border.end(),compare_y_cord);
  for(int i=0;i<border.size();i++){
    temp_pic.set_intensity(border[i].y,border[i].x,0);
    area.push_back(border[i]);
  }
  //temp_pic.print_picture();


  for(int i=0;i<border.size()-1;i++){
    if(temp_pic.get_intensity(border[i].y-1,border[i].x+1) == 0 && border[i].y == border[i+1].y ){
      vector<Point> line;
      line = segment(border[i],border[i+1]);
      for(int j=0;j<line.size();j++){
        temp_pic.set_intensity(line[j].y,line[j].x,0);
        double is_elipse=(double)((line[j].x-center.x)*(line[j].x-center.x)/pow(a,2) + (line[j].y-center.y)*(line[j].y-center.y)/pow(b,2));
        if(is_elipse >= 1){
          area.push_back(Point(line[j].x,line[j].y));
        }
      }
    }
  }
  //temp_pic.print_picture();
  return area;
}

Picture Picture::attenuation_weak_area(vector<Point> area, Point center, unsigned int a, unsigned int b){
  Picture pressure_pic(Mat(y_length,x_length,CV_8UC1,Scalar(255)));
  Picture pressure_buffer(Mat(y_length,x_length,CV_8UC1,Scalar(0)));
  vector<Point> Elipse;
  Elipse=ellipse_nbh(center,a,b);
  for (int i=0;i<x_length;i++){
    for (int j=0;j<y_length;j++){
      float max_distance=max(a,b);
      float r = norm(Point(i*b/max_distance,j*a/max_distance)-Point(center.x*b/max_distance,center.y*a/max_distance));
      float intensity=get_intensity(j,i);
      float transformation=1-exp(-pow((r)/(max_distance),50));
      pressure_buffer.set_intensity(j,i,transformation);
    }
  }

  for (int i=0;i<Elipse.size();i++){
    pressure_pic.set_intensity(Elipse[i].y,Elipse[i].x,get_intensity(Elipse[i].y,Elipse[i].x));
  }



  for (int i=0; i<area.size();i++){
    float intensity=get_intensity(area[i].y,area[i].x);
    float transformation = get_intensity(area[i].y,area[i].x);
    //float transformation=1-(1-intensity)*(1-pressure_buffer.get_intensity(area[i].y,area[i].x));
    pressure_pic.set_intensity(area[i].y,area[i].x,0);
  }
  /*
  for (int i=0;i<Elipse.size();i++){
    float r = norm(Elipse[i]-center);
    float max_distance=max(a,b);
    float intensity=get_intensity(Elipse[i].y,Elipse[i].x);
    float transformation=1-pow((1-intensity),2)*exp(-pow((r+50)/(max_distance),1.5));
    pressure_pic.set_intensity(Elipse[i].y,Elipse[i].x,transformation);

  }
  for (int i=0;i<area.size();i++){
    float r = norm(area[i]-center);
    float max_distance=2*max(a,b);
    float intensity=get_intensity(area[i].y,area[i].x);
    float transformation=1-pow((1-intensity),10)*exp(-pow((r)/(max_distance),4));
    pressure_pic.set_intensity(area[i].y,area[i].x,transformation);

  }*/
  return(pressure_pic);
}
