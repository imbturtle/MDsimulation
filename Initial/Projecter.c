#include<stdio.h>
#include"MD.h"
/********************************************************************************
 Function : use center V for reference.
*********************************************************************************/
double Projector(int N,double **ArrR,double halfL,double R){
    int i,j;
    for(i=0;i<N;i++)for(j=0;j<3;j++)*(*(ArrR+i)+j) = *(*(ArrR+i)+j)*R;
    return halfL*R;
}
