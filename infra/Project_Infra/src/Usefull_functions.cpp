#include "Usefull_functions.h"
#include <iostream>


using namespace std;


float iitof(int color){
  if(color<0 || 255<color){
    cout<<"invalid color"<<"\n";
    return(-1);}
    else return(float(color)/255);

  }

int iftoi(float color){
  if(color<0 || 1<color){
    cout<<"invalid color"<<"\n";
    return(-1);}
    else return(int(color*255));

}

void display_matrix(float ** matrix,int row, int col){
  for (int i=0; i<col;i++){
    for(int j=0; j<row;j++){
      cout << matrix[i][j] << "\t";
    }
    cout << "\n";
  }
}
