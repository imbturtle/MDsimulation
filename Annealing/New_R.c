#include<stdio.h>
#include"Annealing.h"
/****************************************************************************************************************
 Function:R(t+dt) = R(t) + V(t)*dt + 1/6*(4*a(t)+a(t-dt))*dt^2 //beeman O^4
****************************************************************************************************************/
void New_R(int N,double **ArrR,double **ArrV,double **ArrA,double **ArrAP,double **ArrAPP,double halfL,double Dt){ 
	int i, j ;
	for(i=0;i<N;i++)for(j=0;j<3;j++){
		*(*(ArrAPP+i)+j)=*(*(ArrAP+i)+j);
		*(*(ArrAP+i)+j)=*(*(ArrA+i)+j);
		*(*(ArrA+i)+j)=0.0;
	}
	for(i=0;i<N;i++)for(j=0;j<3;j++){
		*(*(ArrR+i)+j) += *(*(ArrV+i)+j) * Dt + 1.0/6.0 * (4.0 * *(*(ArrAP+i)+j) - *(*(ArrAPP+i)+j)) * Dt * Dt ;
        while(*(*(ArrR+i)+j) > halfL) *(*(ArrR+i)+j) -= halfL*2.;
        while(*(*(ArrR+i)+j) <-halfL) *(*(ArrR+i)+j) += halfL*2.;
}}
