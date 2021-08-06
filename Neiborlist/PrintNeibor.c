#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Neiborlist.h"
/***********************************************************************************************************
***********************************************************************************************************/
int PrintNeibor(FILE*fw,int N,int i,int**Nblist,double**ArrR,double**Nbdist,double*halfL){
	int j,l,Nbcount=0;
	double Dij[3];
	extern int Nsearch,Ncut;
	extern double Rsearch,Rcut;
	for(j=0;j<Nsearch;j++){
		if(Nblist[i][j]>=0){
			for(l=0;l<3;l++){
				Dij[l] =  ArrR[Nblist[i][j]][l] - ArrR[i][l];
				while(Dij[l] > halfL[l]) Dij[l] -= halfL[l]*2.0;
				while(Dij[l] <-halfL[l]) Dij[l] += halfL[l]*2.0;
			}
			if(Nbdist[i][j]/Nbdist[i][0]<=Rcut && j<Ncut) fprintf(fw,"%lf %lf %lf\n",Dij[0]/Nbdist[i][0],Dij[1]/Nbdist[i][0],Dij[2]/Nbdist[i][0]);
			if(Nbdist[i][j]/Nbdist[i][0]<=Rcut) Nbcount++;
			if(Nbdist[i][j]/Nbdist[i][0]>Rcut && j<Ncut) fprintf(fw,"%lf %lf %lf\n",0.0,0.0,0.0);
		}
		if(Nblist[i][j]<0 && j<Ncut) fprintf(fw,"%lf %lf %lf\n",0.0,0.0,0.0);
	}
	return Nbcount;
}
