#include<math.h>
#include"Xcate_Ave.h"
/***********************************************************************************************************
***********************************************************************************************************/
double Cal_Nb(int N,int i,int**Nblist,double**ArrR,double**Nbdist,double*halfL){
	int j,k,l;
	double Dij[3],Lij,NbCount=0.0;
	extern int Nearest;
	extern double Rmin;
	for(j=0;j<N;j++){
		if(i==j) continue;
		for(l=0;l<3;l++){
			Dij[l] =  *(*(ArrR+j)+l) - *(*(ArrR+i)+l);
			while(Dij[l] > halfL[l]) Dij[l] -= halfL[l]*2.0;
			while(Dij[l] <-halfL[l]) Dij[l] += halfL[l]*2.0;
		}
		Lij = sqrt(Dij[0] * Dij[0] + Dij[1] * Dij[1] + Dij[2] * Dij[2]);
		if(Lij > Rmin) continue;
		NbCount++;
		for(k=0;k<Nearest;k++){
			if(Lij < Nbdist[i][k]){
				for(l=Nearest;l>k;l--){
					Nbdist[i][l] = Nbdist[i][l-1];
					Nblist[i][l] = Nblist[i][l-1];
				}
				Nbdist[i][k] = Lij;
				Nblist[i][k] = j;
				break;
	}}}
	return NbCount;
}
