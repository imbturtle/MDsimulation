#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"Xcate.h"
double Rmin,qc=0.6;
int Nearest=20,SolidCri=8;
int main(int argc,char*argv[]){
	int i,m,*Connect,**Nblist,matricsize=3,Nparticles,Conf;
	double ***F4,***F6,**Nbdist,**ArrR,**Rq4,**Iq4,**Rq6,**Iq6,*Locq4sq,*Locq6sq,**LocVecq4R,**LocVecq4I,**LocVecq6R,**LocVecq6I\
		,*Locq4,*Locq6,*Locw4,*Locw6,*Si,*Aveq4,*Aveq6,*Avew4,*Avew6,**matric,*eigenvalue,*eigen,**Xcate,halfL[3];
	FILE *fr, *fw;
	clock_t start,end;
	if(argc!=9){
		printf("argv should be 9");
		exit (1);
	}
	fr = fopen(argv[1],"r");
	fw = fopen(argv[2],"w");
	Conf = atoi(argv[3]);
	Nparticles = atoi(argv[4]); 
	Rmin = atof(argv[5]);
	halfL[0] = atof(argv[6]);
	halfL[1] = atof(argv[7]);
	halfL[2] = atof(argv[8]);
	ArrR = Array2nd(Nparticles,3);
	Nbdist = Array2nd(Nparticles,Nearest);
	Nblist = (int **)calloc(Nparticles,sizeof(int *));
	for(i=0;i<Nparticles;i++) Nblist[i] = (int *)calloc(Nearest,sizeof(int));
	F4 = Array3rd(9, 9, 9);
	F6 =  Array3rd(13, 13, 13);
	Rq4 = Array2nd(Nparticles,9);
	Iq4 = Array2nd(Nparticles,9);
	Rq6 = Array2nd(Nparticles,13);
	Iq6 = Array2nd(Nparticles,13);
	Locq4sq = Array1st(Nparticles);
	Locq6sq = Array1st(Nparticles);
	LocVecq4R = Array2nd(Nparticles,9);
	LocVecq4I = Array2nd(Nparticles,9);
	LocVecq6R = Array2nd(Nparticles,13);
	LocVecq6I = Array2nd(Nparticles,13);
	Locq4 = Array1st(Nparticles);
	Locq6 = Array1st(Nparticles);
	Locw4 = Array1st(Nparticles);
	Locw6 = Array1st(Nparticles);
	Si = Array1st(Nparticles);
	Connect = (int *)calloc(Nparticles,sizeof(int));
	Aveq4 = Array1st(Nparticles);
	Aveq6 = Array1st(Nparticles);
	Avew4 = Array1st(Nparticles);
	Avew6 = Array1st(Nparticles);
	matric = Array2nd(4,4);
	eigenvalue = Array1st(4);
	eigen = Array1st(4);
	Xcate = Array2nd(Nparticles,3);
	Winger3j(F4,F6);
	start=clock();
	for(m=0;m<Conf;m++){
		for(i=0;i<Nparticles;i++) fscanf(fr,"%lf %lf %lf",&*(*(ArrR+i)+0),&*(*(ArrR+i)+1),&*(*(ArrR+i)+2));
		Reset(Nparticles,Nblist,Nbdist,Rq4,Iq4,Rq6,Iq6,Locq4sq,Locq6sq,Xcate);
		for(i=0;i<Nparticles;i++) Cal_Nb(Nparticles,i,Nblist,ArrR,Nbdist,halfL);
		for(i=0;i<Nparticles;i++) {
			mat(i,Nblist,halfL,ArrR,matric);
			tred2(matric,matricsize,eigenvalue,eigen);
			tqli(eigenvalue,eigen,matricsize,matric);
			Crystalvec(Nparticles,i,halfL,Nblist,ArrR,LocVecq4R,LocVecq4I,LocVecq6R,LocVecq6I,eigenvalue,matric,Xcate);
			if(i!=0&&i%(Nparticles/5)==0){
				end=clock();
				printf("%d %d cost_time:%ldsec\n",m,i,(end-start)/CLOCKS_PER_SEC);
				start=clock();
		}}
		for(i=0;i<Nparticles;i++) fprintf(fw,"%d %d %lf %lf %lf \n",m,i,Xcate[i][0],Xcate[i][1],Xcate[i][2]);
	}
	return 0;
}



