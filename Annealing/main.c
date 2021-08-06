#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"Annealing.h"
void main(int argc,char* argv[]){
	int i,t,x,CTimes,steps=120000,Nparticles=3000,Balance_steps=100;
	char Ifile[20]="Init_inform.txt";
	double Rcut =2.5,Dt=0.01,Alpha=0.5,T_bath;
	double **ArrR,**ArrV,**ArrA,**ArrAP,**ArrAPP,KE,Potent,AveKE,halfL;
	FILE *frR,*frV,*fwR,*fwV,*fwE,*fwI;
	clock_t start,end;
	if(argc!=8) exit(1) ;
	frR=fopen(argv[1],"r");    
	frV=fopen(argv[2],"r");
	fwR=fopen(argv[3],"w");    
	fwV=fopen(argv[4],"w");
	fwE=fopen(argv[5],"w");
	halfL=atof(argv[6]);
	T_bath=atof(argv[7]);
	fwI=fopen(Ifile,"a");
	///////////////
	ArrR = Array2nd(Nparticles,3);
	ArrV = Array2nd(Nparticles,3);
	for(i=0;i<Nparticles;i++){
		fscanf(frR,"%lf %lf %lf",&*(*(ArrR+i)+0),&*(*(ArrR+i)+1),&*(*(ArrR+i)+2));
		fscanf(frV,"%lf %lf %lf",&*(*(ArrV+i)+0),&*(*(ArrV+i)+1),&*(*(ArrV+i)+2));
	}
	fclose(frR);
	fclose(frV);
	ArrA = Array2nd(Nparticles,3);
	ArrAP = Array2nd(Nparticles,3);
	ArrAPP = Array2nd(Nparticles,3);
	Potent=Acc(Nparticles,ArrR,ArrA,halfL,Rcut);
	KE=Center_V(Nparticles,ArrV);
	///////////
	start=clock();
	for(t=1;t<=steps;t++){
		if(t<=steps/2 && t%Balance_steps==0) Quench(Nparticles,ArrV,KE,T_bath,Alpha);
		New_R(Nparticles,ArrR,ArrV,ArrA,ArrAP,ArrAPP,halfL,Dt);        
    	Potent=Acc(Nparticles,ArrR,ArrA,halfL,Rcut);     
    	New_V(Nparticles,ArrV,ArrA,ArrAP,ArrAPP,Dt);            
		KE=Center_V(Nparticles,ArrV); 
		if(t%(steps/10)==0){
			end=clock();
			printf("T_bath %lf,Temp %lf,Tpotent %lf cost_time:%ldsec\n",T_bath,2.*KE/3./Nparticles,Potent,(end-start)/CLOCKS_PER_SEC);
			start=clock();
		}
		fprintf(fwE,"%d\t%lf\t%lf\t%lf\n",t,KE,Potent,Potent+KE);	
		if((steps-t)<3000)for(i=0;i<Nparticles;i++){
			fprintf(fwR,"%lf %lf %lf\n",*(*(ArrR+i)+0),*(*(ArrR+i)+1),*(*(ArrR+i)+2));
			fprintf(fwV,"%lf %lf %lf\n",*(*(ArrV+i)+0),*(*(ArrV+i)+1),*(*(ArrV+i)+2));
			AveKE += KE;
		}
	}
	if(((2.*AveKE/3./Nparticles/3000)-T_bath)>0.02*T_bath) printf("inbalance");
	fprintf(fwI,"T_bath=%lf\tRecor 3000",T_bath);
}
