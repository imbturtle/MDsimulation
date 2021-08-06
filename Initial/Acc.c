#include<stdio.h>
#include<math.h>
#include"MD.h"
/***************************************************************************
 acc(Rcut, Nparticles, Position bin, box_halflength, Acceleration bin, potential E)
 Function :Lennard Jones Potential 6,12 , V(r)=4*((r)^(-12)+(r)^(-6))+Ar+B
 give acc from potential to particle at dt
 ****************************************************************************/
double Acc(int N,double **ArrR,double **ArrA,double halfL,double rc){
    int i,k,j;
    double r,r2,dis[3]={0.0},r6,r12,potential,Aik,A,B;
    double Potent=0.0;
    A=48.0*pow(rc,-13.0)-24.0*pow(rc,-7.0); //calculate A,B from RC
    B=-52.0*pow(rc,-12.0)+28.0*pow(rc,-6.0);
    for(i=0;i<N-1;i++){
        for(k=i+1;k<N;k++){ //i to k
            r2=0.0 ;
            for(j=0;j<3;j++){
                dis[j]=*(*(ArrR+k)+j)-*(*(ArrR+i)+j) ;
                while(dis[j] > halfL) dis[j] -= halfL*2. ;
                while(dis[j] <-halfL) dis[j] += halfL*2. ;                
                r2 += (dis[j]*dis[j]) ;
            }
            r6=r2*r2*r2 ;
            r12=r6*r6 ;
            r=sqrt(r2) ;
            if(r>rc) potential=0.0 ;
            else potential=4.0*(pow(r,-12.0)-pow(r,-6.0))+A*r+B;
            Potent += potential ;
            for(j=0;j<3;j++){
                if(r>rc) Aik=0.0 ;
                else Aik=24.0*((2.0*dis[j]/(r12*r2))-(dis[j]/(r6*r2)))-A*(dis[j]/r) ; //F=-(delta)V
                *(*(ArrA+i)+j) += -Aik ; //reaction force
                *(*(ArrA+k)+j) += Aik ;  //i to k
	}}}
	return Potent;	
}
