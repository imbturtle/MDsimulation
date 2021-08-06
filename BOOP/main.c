#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"BOOP.h"
int Nearest=20,SolidCri=8;
double Rmin,qc=0.6;
int main(int argc,char*argv[]){
	int i,m,*TanaCrystal,*Connect,**Nblist,Nparticles,Conf;
	double halfL[3],***F4,***F6,**Nbdist,**ArrR,**Rq4,**Iq4,**Rq6,**Iq6,*Locq4sq,*Locq6sq,**LocVecq4R,**LocVecq4I,**LocVecq6R,**LocVecq6I\
		,*Locq4,*Locq6,*Locw4,*Locw6,*Si,*Aveq4,*Aveq6,*Avew4,*Avew6,*NbCount;
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
	NbCount = Array1st(Nparticles);
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
	Connect = (int *)calloc(Nparticles,sizeof(int ));
	Aveq4 = Array1st(Nparticles);
	Aveq6 = Array1st(Nparticles);
	Avew4 = Array1st(Nparticles);
	Avew6 = Array1st(Nparticles);
	TanaCrystal = (int *)calloc(Nparticles,sizeof(int ));
	Winger3j(F4,F6);
	start=clock();
	for(m=0;m<Conf;m++){
		for(i=0;i<Nparticles;i++) fscanf(fr,"%lf %lf %lf ",&*(*(ArrR+i)+0),&*(*(ArrR+i)+1),&*(*(ArrR+i)+2));
		Reset(Nparticles,Nblist,Nbdist,Rq4,Iq4,Rq6,Iq6,Locq4sq,Locq6sq);
		for(i=0;i<Nparticles;i++) NbCount[i]=Cal_Nb(Nparticles,i,Nblist,ArrR,Nbdist,halfL);
		for(i=0;i<Nparticles;i++) Cal_SHP(Nparticles,i,Nblist,ArrR,Rq4,Iq4,Rq6,Iq6,Locq4sq,Locq6sq,LocVecq4R,LocVecq4I,LocVecq6R,LocVecq6I,halfL);
///////////////////////////////Local Bond Orientational Order Parameter/////////////////////////////////////////////
		for(i=0;i<Nparticles;i++) Cal_LocQW(Nparticles,i,Nblist,NbCount,Rq4,Iq4,Rq6,Iq6,F6,F4,Locq4sq,Locq6sq,Locq4,Locq6,Locw4,Locw6);
////////////////////////////cal vec//////////////////////////////////////////////////////////////////////////////////
		for(i=0;i<Nparticles;i++) Si[i]=Cal_SPBOP(Nparticles,i,Connect,Nblist,LocVecq4R,LocVecq4I,LocVecq6R,LocVecq6I);
////////////////////////Average Bond Orientational Order Parameter/////////////////////////////////////////////////////////////
		for(i=0;i<Nparticles;i++) Cal_AveQW(Nparticles,i,Nblist,NbCount,Rq4,Iq4,Rq6,Iq6,F4,F6,Aveq4,Aveq6,Avew4,Avew6);
		for(i=0;i<Nparticles;i++) fprintf(fw,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf\n"\
			,i,Connect[i],Locq4[i],Locq6[i],Locw4[i],Locw6[i],Aveq4[i],Aveq6[i],Avew4[i],Avew6[i]);
		if(m!=0 && m%(Conf/5)==0){
			end=clock();
			printf("%d cost_time:%ldsec\n",m,(end-start)/CLOCKS_PER_SEC);
			start=clock();
	}}
	return 0;
}
