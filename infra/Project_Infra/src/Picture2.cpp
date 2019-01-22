#include "Picture.h"
#include "Usefull_functions.h"
#include<opencv2/opencv.hpp>
#include <math.h>
#include <vector>
#include <algorithm>
using namespace cv;


Picture Picture::cast_rotation_cart(Point centre,float o){//brut method²
  float rad = M_PI/(180)*o;
  Picture res(256,288);
  std::vector<Point> list_res ;
  std::cout << "x_len" << x_length<<'\n';
  std::cout << "y_len" <<y_length<< '\n';
  for(int i = 0  ; i< x_length; i++){
    for(int j = 0; j<y_length; j++){
      int x = int(cos(rad)*(i-centre.x) -sin(rad)*(j-centre.y)) + centre.x;
      int y = int(sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y)) +centre.y;
      // std::cout << "x "<< x << '\n';
      // std::cout << "y" << y<<'\n';
      if((x>0 && x < get_x_len()) && (y>0 && y < get_y_len())){
        res.set_intensity(j,i,get_intensity(y,x));
      }
    }
  }
  return res;


}

std::vector<Point> Picture::original_coordt(float o){
  float rad = M_PI/(180.)*o;
  Point milieu(get_y_len()/2, get_x_len()/2);
  std::vector<Point_<float> > list_res ;
  std::vector<Point> original_coord ;
  for( int i = 0  ; i< milieu.x; i++){
    for( int j = 0; j<milieu.y; j++){
        float x = cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y) + milieu.x;

        float y = sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y) +milieu.y;
        if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
            list_res.push_back(Point2f(x,y));
            original_coord.push_back(Point(i,j));
        }

}

}

  for(int i = milieu.x  ; i< get_x_len(); i++){
    for(int j = 0; j<milieu.y; j++){
      float x = cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y) + milieu.x;
      float y = (sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y))+milieu.y;

      if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
          list_res.push_back(Point2f(x,y));
          original_coord.push_back(Point(i,j));
      }
    }
}
  for(int i = 0 ; i< milieu.x; i++){
    for(int j = milieu.y; j<get_y_len(); j++){
      float x = (cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y)) + milieu.x;
      float y = (sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y))+milieu.y;

      if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
          list_res.push_back(Point2f(x,y));
          original_coord.push_back(Point(i,j));
      }
          }
}
  for(int i = milieu.x  ; i< get_x_len(); i++){
    for(int j = milieu.y; j<get_y_len(); j++){
      float x = cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y)+ milieu.x;
      float y = sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y)+milieu.y;
      if((x < x_length && x >=0 )&& (y < y_length && y >=0) ){
          list_res.push_back(Point2f(x,y));
          original_coord.push_back(Point(i,j));

      }
    }
  }

  return original_coord;}

std::vector<Point_<float>> Picture::list_send_point(float o){
    float rad = M_PI/(180.)*o;
    Point milieu(get_y_len()/2, get_x_len()/2);
    std::vector<Point_<float>> list_res ;
      std::vector<Point_<float>> original_coord ;
    for( int i = 0  ; i< milieu.x; i++){
      for( int j = 0; j<milieu.y; j++){


          float x = cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y) + milieu.x;

          float y = sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y) +milieu.y;
          if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
              list_res.push_back(Point2f(x,y));
              original_coord.push_back(Point(i,j));
          }


}
    }
    for(int i = milieu.x  ; i< get_x_len(); i++){
      for(int j = 0; j<milieu.y; j++){
        float x = cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y) + milieu.x;
        float y = (sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y))+milieu.y;

        if((x < x_length && x>=0 ) && (y < y_length && y >=0) ){
            list_res.push_back(Point2f(x,y));
            original_coord.push_back(Point(i,j));
        }
      }
    }
    for(int i = 0 ; i< milieu.x; i++){
      for(int j = milieu.y; j<get_y_len(); j++){
        float x = (cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y)) + milieu.x;
        float y = (sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y))+milieu.y;
        if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
            list_res.push_back(Point2f(x,y));
            original_coord.push_back(Point(i,j));
        }
            }
}
    for(int i = milieu.x  ; i< get_x_len(); i++){
      for(int j = milieu.y; j<get_y_len(); j++){
        float x = cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y)+ milieu.x;
        float y = sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y)+milieu.y;
        if((x < x_length && x >=0 )&& (y < y_length && y >=0) ){
            list_res.push_back(Point2f(x,y));
            original_coord.push_back(Point(i,j));

        }
      }
    }
    return list_res;}


Picture Picture::bilinear_interpolation_cart( double o){
  Picture res(x_length, y_length); //initialized rotated pic white
  std::vector<Point_<float>> list;
  std::vector<Point> original_coord ;

  float rad = M_PI/(180)*o;
  list = res.list_send_point(o); //list of original pixels
  original_coord = res.original_coordt(o); // list of original coordinates;
//sort if points are in the frame.
  int len = list.size();
  for(int i = 0; i< len; i++){
    if(list[i].x > get_x_len() || list[i].y > get_y_len() || list[i].x<0 || list[i].y <0 ){

      list.erase(list.begin()+i);
      original_coord.erase(original_coord.begin() + i);


      }
  }
  int len_1 = list.size();
        for(int i = 0; i<len; i++){
          Point upright(int(floor(list[i].x)), int(floor(list[i].y)+1));
          Point upleft(int(floor(list[i].x)), int(floor(list[i].y)));
          Point downright(int(floor(list[i].x)+1), int(floor(list[i].y)+1));
          Point downleft(int(floor(list[i].x)+1), int(floor(list[i].y)));
          float intensity_x_0 = ((downleft.x-list[i].x)*get_intensity(upleft.y,upleft.x)) + (list[i].x - upleft.x)*get_intensity(downleft.y,downleft.x);
          float intensity_x_1 = ((downleft.x-list[i].x)*get_intensity(upright.y,upright.x)) + (list[i].x - upleft.x)*get_intensity(downright.y,downright.x);
          float finaaaaaal_intensity = (downright.y-list[i].y)*intensity_x_0 + (list[i].y-upleft.y)*intensity_x_1;
          res.set_intensity(original_coord[i].y,original_coord[i].x,finaaaaaal_intensity );
          }
    return res;

  }


std::vector<Point> Picture::neighb_original_coordt(float o,float a, float b, Point centre){
      float rad = M_PI/(180.)*o;
      Point milieu(get_y_len()/2, get_x_len()/2);
      std::vector<Point_<float> > list_res ;
      std::vector<Point> original_coord ;
      for( int i = 0  ; i< milieu.x; i++){
        for( int j = 0; j<milieu.y; j++){
            if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
            float x = cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y) + milieu.x;

            float y = sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y) +milieu.y;
            if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                list_res.push_back(Point2f(x,y));
                original_coord.push_back(Point(i,j));
            }
          }
    }

    }

      for(int i = milieu.x  ; i< get_x_len(); i++){
        for(int j = 0; j<milieu.y; j++){
          if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
          float x = cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y) + milieu.x;
          float y = (sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y))+milieu.y;

          if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
              list_res.push_back(Point2f(x,y));
              original_coord.push_back(Point(i,j));
          }
        }
      }
    }
      for(int i = 0 ; i< milieu.x; i++){
        for(int j = milieu.y; j<get_y_len(); j++){
            if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
          float x = (cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y)) + milieu.x;
          float y = (sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y))+milieu.y;

          if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
              list_res.push_back(Point2f(x,y));
              original_coord.push_back(Point(i,j));
          }
        }
      }
    }
      for(int i = milieu.x  ; i< get_x_len(); i++){
        for(int j = milieu.y; j<get_y_len(); j++){
            if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
          float x = cos(rad)*(i-milieu.x)- sin(rad)*(j-milieu.y)+ milieu.x;
          float y = sin(rad)*(i-milieu.x)+ cos(rad)*(j-milieu.y)+milieu.y;
          if((x < x_length && x >=0 )&& (y < y_length && y >=0) ){
              list_res.push_back(Point2f(x,y));
              original_coord.push_back(Point(i,j));
            }
          }
        }
      }

      return original_coord;}
std::vector<Point_<float>> Picture::neighb_ellipse_point(float o,float a, float b, Point centre){
            float rad = M_PI/(180.)*o;
            Point milieu(get_y_len()/2, get_x_len()/2);
            std::vector<Point_<float> > list_res ;
            std::vector<Point> original_coord ;
            for( int i = 0  ; i< milieu.x; i++){
              for( int j = 0; j<milieu.y; j++){
                  if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                    float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                    float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;
                  if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                      list_res.push_back(Point2f(x,y));
                      original_coord.push_back(Point(i,j));
                  }
                }
          }

          }

            for(int i = milieu.x  ; i< get_x_len(); i++){
              for(int j = 0; j<milieu.y; j++){
                if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                  float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                  float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;

                if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                    list_res.push_back(Point2f(x,y));
                    original_coord.push_back(Point(i,j));
                }
              }
            }
          }
            for(int i = 0 ; i< milieu.x; i++){
              for(int j = milieu.y; j<get_y_len(); j++){
                  if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                    float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                    float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;

                if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                    list_res.push_back(Point2f(x,y));
                    original_coord.push_back(Point(i,j));
                }
              }
            }
          }
            for(int i = milieu.x  ; i< get_x_len(); i++){
              for(int j = milieu.y; j<get_y_len(); j++){
                  if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                    float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                    float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;
                if((x < x_length && x >=0 )&& (y < y_length && y >=0) ){
                    list_res.push_back(Point2f(x,y));
                    original_coord.push_back(Point(i,j));
                  }
                }
              }
            }

            return list_res;}
Picture Picture::bilinear_interpolation_cart_nghb( double o, Point centre, float a, float b){
              Picture res(x_length, y_length); //initialized rotated pic white
              std::vector<Point_<float>> list;
              std::vector<Point> original_coord ;
              float rad = M_PI/(180)*o;
              list = res.neighb_ellipse_point(-o,a,b,centre); //list of original pixels
              original_coord = res.neighb_original_coordt(-o,a,b,centre); // list of original coordinates;
            //sort if points are in the frame.
              int len = list.size();
              for(int i = 0; i< len; i++){
                if(list[i].x > get_x_len() || list[i].y > get_y_len() || list[i].x<0 || list[i].y <0 ){

                  list.erase(list.begin()+i);
                  original_coord.erase(original_coord.begin() + i);


                  }
              }
              int len_1 = list.size();
                    for(int i = 0; i<len; i++){
                      Point upright(int(floor(list[i].x)), int(floor(list[i].y)+1));
                      Point upleft(int(floor(list[i].x)), int(floor(list[i].y)));
                      Point downright(int(floor(list[i].x)+1), int(floor(list[i].y)+1));
                      Point downleft(int(floor(list[i].x)+1), int(floor(list[i].y)));
                      float intensity_x_0 = ((downleft.x-list[i].x)*get_intensity(upleft.y,upleft.x)) + (list[i].x - upleft.x)*get_intensity(downleft.y,downleft.x);
                      float intensity_x_1 = ((downleft.x-list[i].x)*get_intensity(upright.y,upright.x)) + (list[i].x - upleft.x)*get_intensity(downright.y,downright.x);
                      float finaaaaaal_intensity = (downright.y-list[i].y)*intensity_x_0 + (list[i].y-upleft.y)*intensity_x_1;
                      res.set_intensity(original_coord[i].y,original_coord[i].x,finaaaaaal_intensity );
                      }
                return res;

              }

std::vector<Point_<float>> Picture::swirl_point(float o,float a, float b, Point centre){
                        float rad;
                          Point milieu(get_y_len()/2, get_x_len()/2);
                          std::vector<Point_<float> > list_res ;
                          std::vector<Point> original_coord ;
                          float d_max = sqrt(pow(get_x_len(),2) + pow(get_y_len(),2)) ;
                          for( int i = 0  ; i< milieu.x; i++){
                            for( int j = 0; j<milieu.y; j++){
                                if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                                  float d_ij =sqrt(pow(i-centre.x,2)+pow(j-centre.y,2));
                                  rad = (M_PI*d_ij)/512;
                                  std::cout << "angle =" << rad << '\n';
                                  float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                                  float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;
                                if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                                    list_res.push_back(Point2f(x,y));
                                    original_coord.push_back(Point(i,j));
                                }
                              }
                        }

                        }

                          for(int i = milieu.x  ; i< get_x_len(); i++){
                            for(int j = 0; j<milieu.y; j++){
                              if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                                rad = (M_PI*sqrt(pow(i-centre.x,2)+pow(j-centre.y,2)))/512;
                                float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                                float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;

                              if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                                  list_res.push_back(Point2f(x,y));
                                  original_coord.push_back(Point(i,j));
                              }
                            }
                          }
                        }
                          for(int i = 0 ; i< milieu.x; i++){
                            for(int j = milieu.y; j<get_y_len(); j++){
                                if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                                  rad = (M_PI*sqrt(pow(i-centre.x,2)+pow(j-centre.y,2)))/512;
                                  float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                                  float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;

                              if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                                  list_res.push_back(Point2f(x,y));
                                  original_coord.push_back(Point(i,j));
                              }
                            }
                          }
                        }
                          for(int i = milieu.x  ; i< get_x_len(); i++){
                            for(int j = milieu.y; j<get_y_len(); j++){
                                if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                                  rad = (M_PI*sqrt(pow(i-centre.x,2)+pow(j-centre.y,2)))/512;
                                  float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                                  float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;
                              if((x < x_length && x >=0 )&& (y < y_length && y >=0) ){
                                  list_res.push_back(Point2f(x,y));
                                  original_coord.push_back(Point(i,j));
                                }
                              }
                            }
                          }

                          return list_res;}
std::vector<Point> Picture::original_swirl_point(float o,float a, float b, Point centre){
                                                  float rad = M_PI/(180.)*o;
                                                    Point milieu(get_y_len()/2, get_x_len()/2);
                                                    std::vector<Point_<float> > list_res ;
                                                    std::vector<Point> original_coord ;
                                                    float d_max = sqrt(pow(get_x_len(),2) + pow(get_y_len(),2)) ;
                                                    for( int i = 0  ; i< milieu.x; i++){
                                                      for( int j = 0; j<milieu.y; j++){
                                                          if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                                                            rad = M_PI/(d_max)*sqrt(pow(i-centre.x,2)+pow(j-centre.y,2));
                                                            float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                                                            float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;
                                                          if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                                                              list_res.push_back(Point2f(x,y));
                                                              original_coord.push_back(Point(i,j));
                                                          }
                                                        }
                                                  }

                                                  }

                                                    for(int i = milieu.x  ; i< get_x_len(); i++){
                                                      for(int j = 0; j<milieu.y; j++){
                                                        if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                                                          rad = M_PI/(d_max)*sqrt(pow(i-centre.x,2)+pow(j-centre.y,2));
                                                          float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                                                          float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;

                                                        if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                                                            list_res.push_back(Point2f(x,y));
                                                            original_coord.push_back(Point(i,j));
                                                        }
                                                      }
                                                    }
                                                  }
                                                    for(int i = 0 ; i< milieu.x; i++){
                                                      for(int j = milieu.y; j<get_y_len(); j++){
                                                          if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                                                            rad = M_PI/(d_max)*sqrt(pow(i-centre.x,2)+pow(j-centre.y,2));
                                                            float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                                                            float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;

                                                        if((x < x_length && x>=0 )&& (y < y_length && y >=0) ){
                                                            list_res.push_back(Point2f(x,y));
                                                            original_coord.push_back(Point(i,j));
                                                        }
                                                      }
                                                    }
                                                  }
                                                    for(int i = milieu.x  ; i< get_x_len(); i++){
                                                      for(int j = milieu.y; j<get_y_len(); j++){
                                                          if( (double)((i-centre.x)*(i-centre.x)/pow(a,2) + (j-centre.y)*(j-centre.y)/pow(b,2)) <= 1 ){
                                                            rad = M_PI/(d_max)*sqrt(pow(i-centre.x,2)+pow(j-centre.y,2));
                                                            float x = cos(rad)*(i-centre.x)- sin(rad)*(j-centre.y) + centre.x;
                                                            float y = (sin(rad)*(i-centre.x)+ cos(rad)*(j-centre.y))+centre.y;
                                                        if((x < x_length && x >=0 )&& (y < y_length && y >=0) ){
                                                            list_res.push_back(Point2f(x,y));
                                                            original_coord.push_back(Point(i,j));
                                                          }
                                                        }
                                                      }
                                                    }

                                                    return original_coord;}
Picture Picture::bilinear_interpolation_cart_swirl( double o, Point centre, float a, float b){
                                        Picture res(x_length, y_length); //initialized rotated pic white
                                        std::vector<Point_<float>> list;
                                        std::vector<Point> original_coord ;
                                        float rad;
                                        list = res.swirl_point(-o,a,b,centre); //list of original pixels
                                        original_coord = res.original_swirl_point(-o,a,b,centre); // list of original coordinates;
                                      //sort if points are in the frame.
                                        int len = list.size();
                                        for(int i = 0; i< len; i++){
                                          if(list[i].x > get_x_len() || list[i].y > get_y_len() || list[i].x<0 || list[i].y <0 ){

                                            list.erase(list.begin()+i);
                                            original_coord.erase(original_coord.begin() + i);


                                            }
                                        }
                                        int len_1 = list.size();
                                              for(int i = 0; i<len; i++){
                                                Point upright(int(floor(list[i].x)), int(floor(list[i].y)+1));
                                                Point upleft(int(floor(list[i].x)), int(floor(list[i].y)));
                                                Point downright(int(floor(list[i].x)+1), int(floor(list[i].y)+1));
                                                Point downleft(int(floor(list[i].x)+1), int(floor(list[i].y)));
                                                float intensity_x_0 = ((downleft.x-list[i].x)*get_intensity(upleft.y,upleft.x)) + (list[i].x - upleft.x)*get_intensity(downleft.y,downleft.x);
                                                float intensity_x_1 = ((downleft.x-list[i].x)*get_intensity(upright.y,upright.x)) + (list[i].x - upleft.x)*get_intensity(downright.y,downright.x);
                                                float finaaaaaal_intensity = (downright.y-list[i].y)*intensity_x_0 + (list[i].y-upleft.y)*intensity_x_1;
                                                res.set_intensity(original_coord[i].y,original_coord[i].x,finaaaaaal_intensity );
                                                }
                                          return res;

                                        }

Picture Picture::floating_translation(float o,float tx, float ty){
  Mat image2;
  Mat M=(Mat_<float>(2,3)<<cos(o),-sin(o),tx,sin(o),cos(o),ty);
  warpAffine( picture, image2, M, picture.size());
  return Picture(image2);
}

Picture Picture::erosion_bounded(std::vector<Point> v){
  float min;
  Picture res = clone();
  int len = v.size();
  for(int i = 0; i< len ; i++){
    int x = v[i].x;
    int y = v[i].y;
    std::vector<float> cross;
    for(int j = 0; j< len ; j++){
      if( norm(v[i]-v[j]) <=sqrt(2)){
        cross.push_back(get_intensity(v[j].y,v[j].x));


      }
    }
    min =*std::min_element(cross.begin(), cross.end());
    res.set_intensity(y,x,min);
    cross.clear();
}
    return res;
}


Picture Picture::dilation_bounded(std::vector<Point> v){
  float max;
  Picture res = clone();
  int len = v.size();
  for(int i = 0; i< len ; i++){
    int x = v[i].x;2
    int y = v[i].y;
    std::vector<float> cross;
    for(int j = 0; j< len ; j++){
      if( norm(v[i]-v[j]) <=1){
        cross.push_back(get_intensity(v[j].y,v[j].x));
        std::cout << "intensity" << get_intensity(v[j].y,v[j].x) << '\n' ;
        std::cout << "cross size " << cross.size() <<  '\n';

      }
    }
    max =*std::max_element(cross.begin(), cross.end());
    res.set_intensity(y,x,max);
    cross.clear();
}
    return res;
}

std::vector<Point> Picture::donut(Point centre, int a, int b){
  std::vector<Point> res;
  for(int i = 0; i< get_x_len(); i++){
    for(int j = 0 ; j< get_y_len(); j++){
      Point c(j,i);
      if(norm(c-centre)>a && norm(c-centre)<b){
        res.push_back(c);

      }
    }
  }
return res;
}
