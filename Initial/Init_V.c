#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"MD.h"
/**********************************************************************************************************
 Function : (1)By equal partition theorom 0.5*N*m*v^2 = (3/2)*N*k*T   →   average_v = sqrt(3kT/m)
 average_v = sqrt(T), Vij produce by box muller,   x = sqrt(-2 * log(u)) * cos(2 * M_PI * v) * std + mean;
**********************************************************************************************************/
void Init_V(int N,double **ArrV,double Temp){
	int i,j ;
	double V_0,randA,randB ;
	V_0=sqrt(Temp) ;
	for(i=0;i<N;i++){  //i = Particle#
		for(j=0;j<3;j++){  //j = component
			randA = drand48();
			randB = drand48();
			*(*(ArrV+i)+j) = sqrt(-2 * log(randA)) * cos(2 * M_PI * randB) * V_0;
}}}
