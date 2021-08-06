///////////find BAD && CN/////////////////
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include"BAD.h"
double Rmin,Dtheta=1.0;
int main(int argc, char* argv[]){
	int i,m,Nparticles,Conf,TotalAngle;
	double **ArrR,*Angle;
	double halfL[3],TotalCount=0.0,CNumber[20]={},Neighbor;
	extern double Dtheta;
	clock_t start,end;
	FILE *fr, *fw, *fw1;
	if(argc!=10){
		printf("argv should be 9");
		exit (1);
	}
	fr = fopen(argv[1],"r");
	fw = fopen(argv[2],"w");
	fw1 = fopen(argv[3],"w");
	Conf = atoi(argv[4]);
	Nparticles = atoi(argv[5]);
	Rmin = atof(argv[6]); 
	halfL[0] = atof(argv[7]);
	halfL[1] = atof(argv[8]);
	halfL[2] = atof(argv[9]);
	ArrR = Array2nd(Nparticles,3);
	TotalAngle = (int) 180.0/Dtheta;
	Angle = Array1st(TotalAngle);
	start=clock();
	for(m=0;m<Conf;m++){
		for(i=0;i<Nparticles;i++) fscanf(fr,"%lf %lf %lf",&*(*(ArrR+i)+0),&*(*(ArrR+i)+1),&*(*(ArrR+i)+2));
		for(i=0;i<Nparticles;i++) Cal_BAD(Nparticles,i,halfL,ArrR,Angle,CNumber);	
		if(m!=0 && m%(Conf/5)==0){
			end=clock();
			printf("%d cost_time:%ldsec\n",m,(end-start)/CLOCKS_PER_SEC);
			start=clock();
	}}
	for(i=0;i<TotalAngle;i++) TotalCount += Angle[i];
	for(i=0;i<TotalAngle;i++) fprintf(fw,"%.0f %lf\n",i*Dtheta,Angle[i]/TotalCount);
	for(i=0;i<20;i++) fprintf(fw1,"%d %f\n",i,CNumber[i]/Nparticles/Conf);
	fclose(fr);
	fclose(fw);
	fclose(fw1);
}
