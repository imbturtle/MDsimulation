///////////find BAD && CN/////////////////
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"bad.c"
double Rmin,Dtheta=1.0;
int main(int argc, char* argv[]){
	int i,m,Nparticles,Conf,TotalAngle;
	double **ArrR,*Angle;
	double *L,TotalCount=0.,CNumber[20]={},Neighbor;
	extern double Dtheta;
	FILE *FpR, *FpW, *FpW1;
	if(argc!=7){
		printf("argv should be 5");
		exit (1);
	}
	FpR = fopen(argv[1],"r");
//argv[2] //bcc fcc hcp
	Rmin = atof(argv[3]); 
	FpW = fopen(argv[4],"w");
	FpW1 = fopen(argv[5],"w");
	Conf = atoi(argv[6]);
	Nparticles = ReadN(FpR);
	ArrR = Array2nd(Nparticles,3);
	TotalAngle = (int) 180.0/Dtheta;
	Angle = Array1st(TotalAngle);
	for(m=0;m<Conf;m++){ //m top
		L=ReadData(FpR,ArrR);
		for(i=0;i<Nparticles;i++) CalBAD(Nparticles,i,L,ArrR,Angle,CNumber);	
	}// m tail
	for(i=0;i<TotalAngle;i++) TotalCount += Angle[i];
	for(i=0;i<TotalAngle;i++) fprintf(FpW,"%.0f %lf\n",i*Dtheta,Angle[i]/TotalCount);
	for(i=0;i<20;i++) fprintf(FpW1,"%d %f\n",i,CNumber[i]/Nparticles/Conf);
	fclose(FpR);
	fclose(FpW);
	fclose(FpW1);
}
