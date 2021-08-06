#include<stdio.h>
#include<math.h>
#include<string.h>
#include"BAD.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Cal_BAD(int N,int i,double *halfL,double **ArrR,double *Angle,double *CN){
	int j,k,l,Count=0;
	double Dij[3],Dik[3],Lij,Lik,Rad=180.0/M_PI,TotalCount=0.,NAngle;
	extern double Rmin,Dtheta;
	for(j=0;j<N;j++){
		if(i==j) continue;
		for(l=0;l<3;l++){
			Dij[l] =  ArrR[j][l] - ArrR[i][l];
			while(Dij[l] > halfL[l]) Dij[l] -= halfL[l]*2.0;
			while(Dij[l] <-halfL[l]) Dij[l] += halfL[l]*2.0;
		}
		Lij = sqrt(Dij[0] * Dij[0] + Dij[1] * Dij[1] + Dij[2] * Dij[2]);
		if(Lij > Rmin) continue;
		Count++;
		for(k=j+1;k<N;k++){
			if(i==k) continue;
			for(l=0;l<3;l++){
				Dik[l] =  ArrR[k][l] - ArrR[i][l];
				while(Dik[l] > halfL[l]) Dik[l] -= halfL[l]*2.0;
				while(Dik[l] <-halfL[l]) Dik[l] += halfL[l]*2.0;
			}
			Lik = sqrt(Dik[0] * Dik[0] + Dik[1] * Dik[1] + Dik[2] * Dik[2]);  
			if(Lik > Rmin) continue; 
			NAngle=((acos((Dij[0] * Dik[0] + Dij[1] * Dik[1] + Dij[2] * Dik[2])/(Lij * Lik)) * Rad)/Dtheta);
			if(((Dij[0] * Dik[0] + Dij[1] * Dik[1] + Dij[2] * Dik[2])/(Lij * Lik))<-1.0) NAngle=179.9;
			Angle[(int)NAngle]++;
	}}
	CN[Count]++;
}
