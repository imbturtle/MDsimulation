//**NVT and annealing program with beeman**/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "gaussian.c"

int main(int argc,char* argv[]){
	int i, x, m, j, k, steps=100, Nparticles;
	double **ArrR, **ArrR_0, *L;
	double Density, d0, noise;
	FILE *fprR, *fpwR;
	srand48(time(NULL));
	if(argc != 5){
		printf("argc should be 5.\n");
		exit(1);
	}
	fprR = fopen(argv[1],"r"); //BCC.xyz//FCC.xyz//HCP.xyz
	//argv[2] //bcc fcc hcp
	noise = atof(argv[3]) ;  //0-0.1%
	fpwR = fopen(argv[4],"w");  
	L = Array1st(3);              //Center V bin
	if((strcmp(argv[2],"sc"))==0){
		L[0]=5.0;
		L[1]=5.0;
		L[2]=5.0;
		Nparticles=1000;
	}
	else if((strcmp(argv[2],"bcc"))==0){
		L[0]=6.29961;
		L[1]=6.29961;
		L[2]=6.29961;
		Nparticles=2000;
		d0= 1.0911;
	}
	else if((strcmp(argv[2],"fcc"))==0){
		L[0]=6.3496;
		L[1]=6.3496;
		L[2]=6.3496;
		Nparticles=2048;
		d0= 1.1224;
	}
	else if((strcmp(argv[2],"hcp"))==0){
		L[0]=4.48985;
		L[1]=7.77665;
		L[2]=7.33189;
		Nparticles=2048;
		d0= 1.1224;
	}
	else{
		L[0]=7.2798;
		L[1]=7.2798;
		L[2]=7.2798;
		Nparticles=3000;
	}
/*******************************************************************************/
	ArrR_0 = Array2nd(Nparticles,3);
	ArrR = Array2nd(Nparticles,3);
/*************************initial part*********************************/         
	Density = Nparticles/L[0]/L[1]/L[2];
	printf("%lf %lf\n",L[1],Density);
	for(i=0;i<Nparticles;i++) fscanf(fprR,"%lf %lf %lf",&*(*(ArrR_0+i)+0),&*(*(ArrR_0+i)+1),&*(*(ArrR_0+i)+2));
/****************** ***************************************************************************/
	for(m=0;m<steps;m++){
		New_R( Nparticles, ArrR_0, ArrR, L, noise, d0);
		for(i=0;i<Nparticles;i++)fprintf(fpwR,"%lf %lf %lf\n",*(*(ArrR+i)+0),*(*(ArrR+i)+1),*(*(ArrR+i)+2)); 
		if(m%(steps/5)==0)printf("%d step ok\n",m);
	}
	fclose(fprR);
	fclose(fpwR);
	return 0;
}
