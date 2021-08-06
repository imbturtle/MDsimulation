#include<stdio.h>
#include<math.h>
#include"MD.h"
/********************************************************************************
 Function : Rescale quench
*********************************************************************************/
double Rescale(int N,double **ArrV,double TotalKE,double Target){
    int i,j;
    double F;
    F = sqrt(Target/TotalKE);
    for(i=0;i<N;i++)for(j=0;j<3;j++)*(*(ArrV+i)+j) = *(*(ArrV+i)+j)*F;
    return 0.;
}
