#include<stdlib.h>
#include<math.h>
#include"Neiborlist.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Cal_Nb(int N,int i,int**Nblist,double**ArrR,double**Nbdist,double*halfL){
	int j,k,l;
	double Dij[3],Lij;
	extern double Rsearch;
	extern int Nsearch;
	for(j=0;j<N;j++){
		if(i==j) continue;
		for(l=0;l<3;l++){
			Dij[l] =  *(*(ArrR+j)+l) - *(*(ArrR+i)+l);
			while(Dij[l] > halfL[l]) Dij[l] -= halfL[l]*2.0;
			while(Dij[l] <-halfL[l]) Dij[l] += halfL[l]*2.0;
		}
		Lij = sqrt(Dij[0] * Dij[0] + Dij[1] * Dij[1] + Dij[2] * Dij[2]);
		if(Lij > Rsearch) continue;
		for(k=0;k<Nsearch;k++)if(Lij < Nbdist[i][k]){
			for(l=Nsearch;l>k;l--){
				Nbdist[i][l] = Nbdist[i][l-1];
				Nblist[i][l] = Nblist[i][l-1];
			}
			Nbdist[i][k] = Lij;
			Nblist[i][k] = j;
			break;
}}}
