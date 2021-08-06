#include<stdio.h>
#include<math.h>
#include"Annealing.h"
/********************************************************************************
 Function : use center V for reference.
*********************************************************************************/
void Quench(int N,double **ArrV,double KE,double Tb,double Alpha){
	int i,j;
	double ksi;
	ksi=(1.0+Alpha*(sqrt(Tb/(2.0*KE/3.0/N))-1.0));
	for(i=0;i<N;i++)for(j=0;j<3;j++) *(*(ArrV+i)+j) *= ksi;
}
