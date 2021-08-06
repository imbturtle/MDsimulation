#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include"RDF.h"
double dr=0.01;
int main(int argc,char* argv[]){
	int i,m,Nparticles,Conf,Totaldr,Halfdr;
	double **ArrR,*g,*Gofr,halfL[3];
	FILE *fr, *fw;
	clock_t start,end;
	if(argc!=8){
		printf("argv should be 8");
		exit (1);
	}
	fr = fopen(argv[1],"r");
	fw = fopen(argv[2],"w");
	Conf = atoi(argv[3]);
	Nparticles = atoi(argv[4]);
	halfL[0] = atof(argv[5]);
	halfL[1] = atof(argv[6]);
	halfL[2] = atof(argv[7]); 
	Totaldr = (int)(30./dr);	
	Halfdr = (int)(Totaldr/2.0);
	g = Array1st(Totaldr);
	Gofr = Array1st(Totaldr);
	ArrR = Array2nd(Nparticles,3);
	start=clock();
	for(m=0;m<Conf;m++){
		for(i=0;i<Nparticles;i++) fscanf(fr,"%lf %lf %lf ",&*(*(ArrR+i)+0),&*(*(ArrR+i)+1),&*(*(ArrR+i)+2));
		for(i=0;i<Nparticles;i++) Cal_NBD(Nparticles,i,halfL,ArrR,g);
		Cal_RDF(Nparticles,Halfdr,halfL,g,Gofr);
		if(m!=0 && m%(Conf/5)==0){
			end=clock();
			printf("%d cost_time:%ldsec\n",m,(end-start)/CLOCKS_PER_SEC);
			start=clock();
	}}
	for(i=0;i<Halfdr;i++) fprintf(fw,"%.2lf %.2lf\n",(dr/2.0)+(double)i*dr,*(Gofr+i)/Nparticles/Conf);
	fclose(fr);
	fclose(fw);
	return 0;
}
