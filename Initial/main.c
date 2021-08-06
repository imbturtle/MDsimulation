#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"MD.h"

void main(){
	int i,t,x,CTimes,steps=1000,Nparticles=3000,Balance_steps=1;
	char Rfile[20]="Rinit.txt",Vfile[20]="Vinit.txt",Ifile[20]="Init_inform.txt";
	double Rcut =2.5,Dt=0.01,Density_tar=0.972,Density_init=0.5,Compress_rate=0.97,Temp=1.0;
	double **ArrR,**ArrV,**ArrA,**ArrAP,**ArrAPP,TotalKE=0,KE,Potent,halfL,Target_KE,A,B;
	FILE *fwR,*fwV,*fwI;
	clock_t start,end;
	fwR=fopen(Rfile,"w");    
	fwV=fopen(Vfile,"w");	
	fwI=fopen(Ifile,"w");
	ArrR = Array2nd(Nparticles,3);
	ArrV = Array2nd(Nparticles,3);
	ArrA = Array2nd(Nparticles,3);
	ArrAP = Array2nd(Nparticles,3);
	ArrAPP = Array2nd(Nparticles,3);
	halfL=BoxLeng(Nparticles,Density_init);
	CTimes = (int)((log(Density_init/Density_tar))/log(Compress_rate))/3.0;
	Target_KE = 3.0*Nparticles*0.5*Temp*Balance_steps;
	Init_R(Nparticles,ArrR,halfL);
	Init_V(Nparticles,ArrV,Temp);
	Potent=Acc(Nparticles,ArrR,ArrA,halfL,Rcut);
	start=clock();
	for(x=0;x<=CTimes+1;x++){
		for(t=1;t<=steps;t++){
			New_R(Nparticles,ArrR,ArrV,ArrA,ArrAP,ArrAPP,halfL,Dt);     
    	    Potent=Acc(Nparticles,ArrR,ArrA,halfL,Rcut);    
    	    New_V(Nparticles,ArrV,ArrA,ArrAP,ArrAPP,Dt);            
			KE=Center_V(Nparticles,ArrV);
			TotalKE+=KE;
			if(t%(steps/5)==0){
				end=clock();
				printf("%d %d %d Tem %lf,Potent %lf cost_time:%ldsec\n",CTimes,x,t,2.*KE/3./Nparticles,Potent,(end-start)/CLOCKS_PER_SEC);
				start=clock();
			}
			if(x<=CTimes && t%Balance_steps==0) TotalKE=Rescale(Nparticles,ArrV,TotalKE,Target_KE);
		}
		if(x<CTimes)halfL=Projector(Nparticles,ArrR,halfL,Compress_rate);
		if(x==CTimes)halfL=Projector_F(Nparticles,ArrR,halfL,Density_tar);
	}
	for(i=0;i<Nparticles;i++){
		fprintf(fwR,"%lf %lf %lf\n",*(*(ArrR+i)+0),*(*(ArrR+i)+1),*(*(ArrR+i)+2));
		fprintf(fwV,"%lf %lf %lf\n",*(*(ArrV+i)+0),*(*(ArrV+i)+1),*(*(ArrV+i)+2));
	}
	A = 48.0*pow(Rcut,-13.0)-24.0*pow(Rcut,-7.0);
	B =-52.0*pow(Rcut,-12.0)+28.0*pow(Rcut,-6.0);
	fprintf(fwI,"N=%d\nhalfL=%lf\nRcut=%lf\nLJ fixed term=\nA:%lf\nB:%lf",Nparticles,halfL,Rcut,A,B);
}
