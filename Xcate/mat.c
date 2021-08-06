#include"Xcate.h"
/***********************************************************************************************************
***********************************************************************************************************/
void mat(int i,int**Nblist,double*halfL,double**ArrR,double**matric){
	int o,j,l;
	double Dij[3];
	extern int Nearest;
	for(o=0;o<4;o++)for(l=0;l<4;l++) matric[o][l]=0.0;	
	for(j=0;j<Nearest;j++){
		if(Nblist[i][j]<0) break;
		for(l=0;l<3;l++){
			Dij[l] =  ArrR[Nblist[i][j]][l] - ArrR[i][l];
			while(Dij[l] > halfL[l]) Dij[l] -= halfL[l]*2.0;
			while(Dij[l] <-halfL[l]) Dij[l] += halfL[l]*2.0;
		}
		matric[1][1]+=(Dij[1]*Dij[1]+Dij[2]*Dij[2]);
		matric[2][2]+=(Dij[0]*Dij[0]+Dij[2]*Dij[2]);
		matric[3][3]+=(Dij[0]*Dij[0]+Dij[1]*Dij[1]);
		matric[1][2]+=(-Dij[0]*Dij[1]);
		matric[1][3]+=(-Dij[0]*Dij[2]);
		matric[2][3]+=(-Dij[1]*Dij[2]);
	}
	matric[2][1]=matric[1][2];
	matric[3][1]=matric[1][3];
	matric[3][2]=matric[2][3];
}
