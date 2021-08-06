#include<math.h>
#include"Xcate.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Cal_CryNb(int N,int i,int**Nblist,double**ArrR,double*halfL,double Rcut){
	int j,k=0,l;
	double Dij[3],Lij;
	for(j=0;j<N;j++){
		if(i==j) continue;
		for(l=0;l<3;l++){
			Dij[l] =  *(*(ArrR+j)+l) - *(*(ArrR+i)+l);
			while(Dij[l] > halfL[l]) Dij[l] -= halfL[l]*2.0;
			while(Dij[l] <-halfL[l]) Dij[l] += halfL[l]*2.0;
		}
		Lij = sqrt(Dij[0] * Dij[0] + Dij[1] * Dij[1] + Dij[2] * Dij[2]);
		if(Lij > Rcut) continue;
		Nblist[i][k]=j;
		k++;
}}

