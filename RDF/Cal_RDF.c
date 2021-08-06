#include<stdio.h>
#include<math.h>
#include<string.h>
#include"RDF.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Cal_RDF(int N,double Halfdr,double*halfL,double*g,double*G){
	int i;
	double Dens;
	extern double dr;
	Dens=N/halfL[0]/halfL[1]/halfL[2]/8.0;
	for(i=0;i<Halfdr;i++) *(g+i) /= (4.0/3.0)*M_PI*(pow(((i+1)*dr),3.0)-pow(i*dr,3.0))*Dens;
	for(i=0;i<Halfdr;i++) *(G+i) += *(g+i);
	memset(g,0,sizeof(double)*Halfdr);
}

