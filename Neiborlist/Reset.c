#include<stdlib.h>
#include<string.h>
#include"Neiborlist.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Reset(int N,int **Nblist,double **Nbdist){
	int i,j;
	extern int Nsearch;
	for(i=0;i<N;i++)for(j=0;j<Nsearch;j++){
		Nbdist[i][j] = 99.0;
		Nblist[i][j] = -1;
}}
