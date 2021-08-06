#include<stdlib.h>
#include<math.h>
#include"BOOP.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Cal_LocQW(int N,int i,int**Nblist,double*NbCount,double**Rq4,double**Iq4,double**Rq6,double**Iq6,double***F6,double***F4\
	,double*Locq4sq,double*Locq6sq,double*Locq4,double*Locq6,double*Locw4,double*Locw6){
	int j,m,m1,m2,m3;
	double w4top=0.0,w6top=0.0;
	for(m1=0;m1<13;m1++)for(m2=0;m2<13;m2++)for(m3=0;m3<13;m3++){
		if(m1+m2+m3==12 && m1<9 && m2<9 && m3<9){
			w4top += F4[m1][m2][m3] * (Rq4[i][m1]*Rq4[i][m2]*Rq4[i][m3] - Iq4[i][m1]*Iq4[i][m2]*Rq4[i][m3] -\
				Rq4[i][m1]*Iq4[i][m2]*Iq4[i][m3] - Iq4[i][m1]*Rq4[i][m2]*Iq4[i][m3]);//Real part
			/*w4topi += F4[m1][m2][m3] * (-Rq4[i][m1]*Rq4[i][m2]*Iq4[i][m3] + Iq4[i][m1]*Iq4[i][m2]*Iq4[i][m3] +\
				Rq4[i][m1]*Iq4[i][m2]*Rq4[i][m3] + Iq4[i][m1]*Rq4[i][m2]*Rq4[i][m3]);//image part*/
		}
		if(m1+m2+m3==18){
			w6top += F6[m1][m2][m3] * (Rq6[i][m1]*Rq6[i][m2]*Rq6[i][m3] - Iq6[i][m1]*Iq6[i][m2]*Rq6[i][m3] -\
				Rq6[i][m1]*Iq6[i][m2]*Iq6[i][m3] - Iq6[i][m1]*Rq6[i][m2]*Iq6[i][m3]);//Real part
			/*w6topi += F6[m1][m2][m3] * (-Rq6[i][m1]*Rq6[i][m2]*Iq6[i][m3] + Iq6[i][m1]*Iq6[i][m2]*Iq6[i][m3] +\
				Rq6[i][m1]*Iq6[i][m2]*Rq6[i][m3] + Iq6[i][m1]*Rq6[i][m2]*Rq6[i][m3]);//image part*/
	}}
	Locq4[i] = sqrt((4.0*M_PI)/(9.0)*Locq4sq[i])/NbCount[i];
	Locq6[i] = sqrt((4.0*M_PI)/(13.0)*Locq6sq[i])/NbCount[i];
	Locw4[i] = w4top/pow(Locq4sq[i],1.5);
	Locw6[i] = w6top/pow(Locq6sq[i],1.5);
}
