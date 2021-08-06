#include<stdio.h>
#include"MD.h"
/************************************************************************************************************
 Function:check the particle is overlapping or not (here dis<1 is not sensible)
 Return: 1 (dis<1 redo init_R) or 0 ()
*************************************************************************************************************/
int Chek_init_R(double **ArrR,double halfL,int I){
	int j,m;
	double d,dis;
	for(m=0;m<I;m++){
		d=0.0 ;
		for(j=0;j<3;j++){
			dis = *(*(ArrR+I)+j) - *(*(ArrR+m)+j) ; 
			while(dis > halfL) dis -= halfL*2.;
			while(dis <-halfL) dis += halfL*2.;
			d += (dis*dis) ;
		}
	if(d<1.0)return 1 ;
	}
	return 0 ;
}
