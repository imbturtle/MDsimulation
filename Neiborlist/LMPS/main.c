#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include"neiborlist.c"
double Rsearch,Rcut;
int Nsearch,Ncut;
int main(int argc,char* argv[]){
	int i,m,Nparticles,Conf,**Nblist,NbCount;
	double **ArrR,**Nbdist,*L,Nbsta[900]={};
	FILE *FpR, *FpW1, *FpW2;
	if(argc != 8){
		printf("argc should be 8.\n");
		exit(1);
	}
	FpR = fopen(argv[1],"r");
	FpW1 = fopen(argv[3],"w"); 
	FpW2 = fopen(argv[4],"w");
	Conf = atoi(argv[5]);
	Ncut = atoi(argv[6]);
	Rcut = atof(argv[7]);
	Rsearch = Rcut*3.0;
	Nsearch = Ncut*3;
	Nparticles = ReadNLAM(FpR);
/*******************************************************************************/
	ArrR = Array2nd(Nparticles,3);
	Nbdist = Array2nd(Nparticles,Nsearch);
	Nblist = (int **)calloc(Nparticles,sizeof(int *));
	for(i=0;i<Nparticles;i++) Nblist[i] = (int *)calloc(Nsearch,sizeof(int));
/****************** ***************************************************************************/
	for(m=0;m<Conf;m++){
		L=ReadDataLAM(FpR,ArrR);
		Reset(Nparticles,Nblist,Nbdist);
		for(i=0;i<Nparticles;i++) {
			CalNb(Nparticles,i,Nblist,ArrR,Nbdist,L);
			NbCount = PrintNeibor(FpW1,Nparticles,i,Nblist,ArrR,Nbdist,L);
			Nbsta[NbCount]++;
		}
		if(m>20)if(m%20==0) printf("m %d ok",m);		
	}
	for(i=0;i<200;i++)fprintf(FpW2,"%d %.2lf\n",i,Nbsta[i]/Nparticles/Conf);
/************************************************************************************************/
	fclose(FpR);
	fclose(FpW1);
	fclose(FpW2);
	return 0;
}
