#include<string.h>
#include"BOOP.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Reset(int N,int**Nblist,double**Nbdist,double**Rq4,double**Iq4,double**Rq6,double**Iq6,double*Locq4sq,double*Locq6sq){
	int i,j,m;
	for(i=0;i<N;i++)for(j=0;j<12;j++){
		Nbdist[i][j] = 99.0;
		Nblist[i][j] = -1;
		for(m=0;m<13;m++){
			if(m<9){
				Rq4[i][m] = 0.0;
				Iq4[i][m] = 0.0;
			}
			Rq6[i][m] = 0.0;
			Iq6[i][m] = 0.0;
	}}
	memset(Locq4sq,0, sizeof(double)*N);
	memset(Locq6sq,0, sizeof(double)*N);
}
