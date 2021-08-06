#include<stdio.h>
#include<math.h>
#include"MD.h"
/********************************************************************************
 Function : use center V for reference.
*********************************************************************************/
double Projector_F(int N,double **ArrR,double halfL,double Target){
    int i,j;
    double R;
    R = pow((N/Target),1.0/3.0)/halfL/2.0;
    for(i=0;i<N;i++)for(j=0;j<3;j++)*(*(ArrR+i)+j) = *(*(ArrR+i)+j)*R;
    return halfL*R;
}
