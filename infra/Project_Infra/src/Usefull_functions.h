#ifndef USEFULL_FUNCTIONS_H
#define USEFULL_FUNCTIONS_H





/*  ------
iitof : Takes color intensity as an int between [0,255] and converts it into a float between [0,1]
----     */
float iitof(int color);

/*  ------
iftoi : Takes color intensity a float between [0,255] and converts it into an int between [0,1]
----     */
int iftoi(float color);

/*  ------
display_mat : Takes a matrix, it's number of row, it's number of col and display it
----     */
void display_matrix(float ** matrix,int row, int col);

























#endif
