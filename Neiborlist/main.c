#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include"Neiborlist.h"
double Rsearch,Rcut;
int Nsearch,Ncut;
int main(int argc,char* argv[]){
	int i,m,Nparticles,Conf,**Nblist,NbCount;
	double **ArrR,**Nbdist,halfL[3],Nbsta[900]={};
	FILE *fr, *fw, *fw1;
	clock_t start,end;
	if(argc != 11){
		printf("argc should be 11.\n");
		exit(1);
	}
	fr = fopen(argv[1],"r");
	fw = fopen(argv[2],"w"); 
	fw1 = fopen(argv[3],"w");
	Conf = atoi(argv[4]);
	Nparticles = atoi(argv[5]);
	halfL[0] = atof(argv[6]);
	halfL[1] = atof(argv[7]);
	halfL[2] = atof(argv[8]); 
	Ncut = atoi(argv[9]);
	Rcut = atof(argv[10]);
	Rsearch = Rcut*3.0;
	Nsearch = Ncut*3;
/*******************************************************************************/
	ArrR = Array2nd(Nparticles,3);
	Nbdist = Array2nd(Nparticles,Nsearch);
	Nblist = (int **)calloc(Nparticles,sizeof(int *));
	for(i=0;i<Nparticles;i++) Nblist[i] = (int *)calloc(Nsearch,sizeof(int));
/****************** ***************************************************************************/
	start=clock();
	for(m=0;m<Conf;m++){
		for(i=0;i<Nparticles;i++) fscanf(fr,"%lf %lf %lf ",&*(*(ArrR+i)+0),&*(*(ArrR+i)+1),&*(*(ArrR+i)+2));
		Reset(Nparticles,Nblist,Nbdist);
		for(i=0;i<Nparticles;i++) {
			Cal_Nb(Nparticles,i,Nblist,ArrR,Nbdist,halfL);
			NbCount = PrintNeibor(fw,Nparticles,i,Nblist,ArrR,Nbdist,halfL);
			Nbsta[NbCount]++;
		}
		if(m!=0 && m%(Conf/5)==0){
			end=clock();
			printf("%d cost_time:%ldsec\n",m,(end-start)/CLOCKS_PER_SEC);
			start=clock();
	}}
	for(i=0;i<200;i++)fprintf(fw1,"%d %.2lf\n",i,Nbsta[i]/Nparticles/Conf);
/************************************************************************************************/
	fclose(fr);
	fclose(fw);
	fclose(fw1);
	return 0;
}
