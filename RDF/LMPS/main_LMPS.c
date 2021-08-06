#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"RDF_LMPS.h"
double dr=0.01;
int main(int argc,char* argv[]){
	int i,m,Nparticles,Conf,Totaldr,Halfdr;
	double **ArrR,*g,*Gofr,*halfL;
	FILE *fr, *fw;
	if(argc!=4){
		printf("argv should be 4");
		exit (1);
	}
	fr = fopen(argv[1],"r");
	fw = fopen(argv[2],"w");
	Conf = atoi(argv[3]); 
	Totaldr = (int)(30./dr);	
	Halfdr = (int)(Totaldr/2.0);
	g = Array1st(Totaldr);
	Gofr = Array1st(Totaldr);
	Nparticles = Read_N(fr);
	ArrR = Array2nd(Nparticles,3);
	for(m=0;m<Conf;m++){//m top
		halfL=Read_Data(fr,ArrR);
		for(i=0;i<Nparticles;i++) Cal_NBD(Nparticles,i,halfL,ArrR,g);
		Cal_RDF(Nparticles,Halfdr,halfL,g,Gofr);
	}//m tail
	for(i=0;i<Halfdr;i++) fprintf(fw,"%.3lf %.3lf\n",(dr/2.0)+(double)i*dr,*(Gofr+i)/Nparticles/Conf);
	fclose(fr);
	fclose(fw);
	return 0;
}
