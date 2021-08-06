#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"RDF_LMPS.h"
/***********************************************************************************************************
Read lmps conf R L
***********************************************************************************************************/
double *Read_Data(FILE*fr,double**ArrR){
	char garba[100];
	int i,N;
	double RI[3],RF[3];
	static double halfL[3];
	fscanf(fr,"%[^\n] ",garba);
	fscanf(fr,"%[^\n] ",garba);
	fscanf(fr,"%[^\n] ",garba);
	fscanf(fr,"%[^\n] ",garba);
	N = atoi(garba);
	fscanf(fr,"%[^\n] ",garba);
	fscanf(fr,"%[^ ] ",garba);
	RI[0] = atof(garba);
	fscanf(fr,"%[^\n] ",garba);
	RF[0] = atof(garba);
	fscanf(fr,"%[^ ] ",garba);
	RI[1] = atof(garba);
	fscanf(fr,"%[^\n] ",garba);
	RF[1] = atof(garba);
	fscanf(fr,"%[^ ] ",garba);
	RI[2] = atof(garba);
	fscanf(fr,"%[^\n] ",garba);
	RF[2] = atof(garba);	
	fscanf(fr,"%[^\n] ",garba);
	for(i=0;i<3;i++) halfL[i] = fabs(RF[i]-RI[i])/2.;
	for(i=0;i<N;i++) fscanf(fr,"%lf %lf %lf ",&*(*(ArrR+i)+0),&*(*(ArrR+i)+1),&*(*(ArrR+i)+2));
	return halfL;
}
