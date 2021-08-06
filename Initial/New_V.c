#include<stdio.h>
#include"MD.h"
/*******************************************************************************************************************
 Function:V(t+dt)*dt = V(t) + 1/6*(2*a(t+dt)+5a(t)-a(t-dt)))*dt  beeman predict O^3
********************************************************************************************************************/
void New_V(int N,double **ArrV,double **ArrA,double **ArrAP,double **ArrAPP,double Dt){
    int i ,j ;
    for(i=0;i<N;i++)for(j=0;j<3;j++){
    	*(*(ArrV+i)+j) += 1.0/6.0 * (2.0 * *(*(ArrA+i)+j) + 5.0 * *(*(ArrAP+i)+j) - *(*(ArrAPP+i)+j)) * Dt ; 
}}
