#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"MD.h"
/************************************************************************************************************
 Funtion : Random value within [-L/2 , L/2]
 Box muller
************************************************************************************************************/
void Init_R(int N,double **ArrR,double halfL){
	int i,j;
	srand48(1) ;//srand48(seed) random value for drand,seed can be int or time(NULL)
	for(i=0;i<N;i++){
		do{
			for(j=0;j<3;j++)*(*(ArrR+i)+j) = halfL * (1.0-2.0*drand48()) ;
		}
		while(Chek_init_R(ArrR,halfL,i)) ; /*dis<1 is not sencible(true non-zero) re-do*/
}}
