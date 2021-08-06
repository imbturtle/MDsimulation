#include<stdio.h>
#include<stdlib.h>
#include"MD.h"
/**********************************************************************************************************
Array def
**********************************************************************************************************/

double *Array1st(int col){
	double *Name;
	Name = (double *)calloc(col,sizeof(double));
	return Name;
}
double **Array2nd(int col, int row){
	int i;
	double **Name;
	Name = (double **)malloc(col*sizeof(double *));
	for(i=0;i<col;i++){
		Name[i] = (double *)malloc(row*sizeof(double));
	}
	return Name;
}
double ***Array3rd(int col, int row,int hight){
	int i,j;
	double ***Name;
	Name = (double ***)malloc(col*sizeof(double **));
	for(i=0;i<col;i++){
		Name[i] = (double **)malloc(row*sizeof(double *));
		for(j=0;j<row;j++) Name[i][j] = (double *)calloc(hight,sizeof(double));
	}
	return Name;
}
