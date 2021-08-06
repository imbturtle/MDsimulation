#include<stdio.h>
#include<stdlib.h>
#include"RDF_LMPS.h"
/***********************************************************************************************************
Read lmps conf N
***********************************************************************************************************/
int Read_N(FILE*fr){
	char garba[100];
	int N;
	fscanf(fr,"%[^\n] ",garba);
	fscanf(fr,"%[^\n] ",garba);
	fscanf(fr,"%[^\n] ",garba);
	fscanf(fr,"%[^\n] ",garba);
	N = atoi(garba);
	rewind(fr);
	return N;
}
