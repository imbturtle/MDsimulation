#include<stdio.h>
#include<math.h>
#include"RDF_LMPS.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Cal_NBD(int N,int i,double*halfL,double**ArrR,double*g){
	int j,k,n;
	double Dij[3];
	extern double dr;
	for(j=0;j<N;j++){
		if(i==j) continue;
		for(k=0;k<3;k++){
			Dij[k] = *(*(ArrR+j)+k) - *(*(ArrR+i)+k);
			while(Dij[k] >  halfL[k]) Dij[k] -= halfL[k]*2.0;
			while(Dij[k] < -halfL[k]) Dij[k] += halfL[k]*2.0;
		}
		n = (int)(sqrt(Dij[0]*Dij[0]+Dij[1]*Dij[1]+Dij[2]*Dij[2])/dr);  
		g[n]++;
}}

