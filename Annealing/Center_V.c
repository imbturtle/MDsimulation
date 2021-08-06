#include<stdio.h>
#include"Annealing.h"
/********************************************************************************
 Function : use center V for reference.
*********************************************************************************/
double Center_V(int N,double **ArrV){
    int i,j;
    double CV[3]={0.0},KE;
    for(j=0;j<3;j++){
        for(i=0;i<N;i++) *(CV+j) += *(*(ArrV+i)+j);
        *(CV+j) /= N;
    }
    for(i=0;i<N;i++)for(j=0;j<3;j++){
    	*(*(ArrV+i)+j) -= CV[j] ;
    	KE += 0.5 * (*(*(ArrV+i)+j) * *(*(ArrV+i)+j));
	}
	return KE;
}
