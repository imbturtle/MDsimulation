#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include"boo.c"
#define pi M_PI
double Rmin;
int main(int argc,char*argv[]){
	int i,m,*TanaCrystal,*Connect,**Nblist,Nparticles,Conf;//parameter
	double *L,***F4,***F6,**Nbdist,**ArrR,**Rq4,**Iq4,**Rq6,**Iq6,*Locq4sq,*Locq6sq,**LocVecq4R,**LocVecq4I,**LocVecq6R,**LocVecq6I\
		,*Locq4,*Locq6,*Locw4,*Locw6,*Si,*Aveq4,*Aveq6,*Avew4,*Avew6,*NbCount;
	FILE *FpR, *FpW1;
	if(argc!=6){
		printf("argv should be 5");
		exit (1);
	}
	FpR = fopen(argv[1],"r"); 
//argv[2] //bcc fcc hcp
	Rmin = atof(argv[3]);
	FpW1 = fopen(argv[4],"w"); 
	Conf = atoi(argv[5]); 
	Nparticles = ReadN(FpR);
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
	for(m=0;m<Conf;m++){
		L=ReadData(FpR,ArrR);
		Reset(Nparticles,Nblist,Nbdist,Rq4,Iq4,Rq6,Iq6,Locq4sq,Locq6sq);
		for(i=0;i<Nparticles;i++) NbCount[i]=CalNb(Nparticles,i,Nblist,ArrR,Nbdist,L);
		for(i=0;i<Nparticles;i++) CalSHP(Nparticles,i,Nblist,ArrR,Rq4,Iq4,Rq6,Iq6,Locq4sq,Locq6sq,LocVecq4R,LocVecq4I,LocVecq6R,LocVecq6I,L);
///////////////////////////////Local Bond Orientational Order Parameter/////////////////////////////////////////////
		for(i=0;i<Nparticles;i++) CalLocQW(Nparticles,i,Nblist,NbCount,Rq4,Iq4,Rq6,Iq6,F6,F4,Locq4sq,Locq6sq,Locq4,Locq6,Locw4,Locw6);
////////////////////////////cal vec//////////////////////////////////////////////////////////////////////////////////
		for(i=0;i<Nparticles;i++) Si[i]=CalSPBOP(Nparticles,i,Connect,Nblist,LocVecq4R,LocVecq4I,LocVecq6R,LocVecq6I);
//		for(i=0;i<Nparticles;i++)fprintf(FpW,"%d %d %lf\n",m,i,Si[i]);
////////////////////////Average Bond Orientational Order Parameter/////////////////////////////////////////////////////////////
		for(i=0;i<Nparticles;i++)CalAveQW(Nparticles,i,Nblist,NbCount,Rq4,Iq4,Rq6,Iq6,F4,F6,Aveq4,Aveq6,Avew4,Avew6);
		for(i=0;i<Nparticles;i++)fprintf(FpW1,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf\n"\
			,i,Connect[i],Locq4[i],Locq6[i],Locw4[i],Locw6[i],Aveq4[i],Aveq6[i],Avew4[i],Avew6[i]);
//		if(m%(Conf/5)==0)printf("%d conf ok",m);
	}// m tail
	return 0;
}



