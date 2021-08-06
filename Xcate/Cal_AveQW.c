#include<math.h>
#include"Xcate.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Cal_AveQW(int N,int i,int**Nblist,double*NbCount,double**Rq4,double**Iq4,double**Rq6,double**Iq6\
	,double***F4,double***F6,double*Aveq4,double*Aveq6,double*Avew4,double*Avew6){
	int j,k,m,m1,m2,m3;
	double Aveq4mR[9],Aveq4mI[9],Aveq6mR[13],Aveq6mI[13],Aveq4sq=0.0,Aveq6sq=0.0,\
		w4top=0.0,w6top=0.0,Count=0.0,w6topi,w4topi;
	extern int Nearest;
	for(m=0;m<13;m++){
		if(m<9){
			Aveq4mR[m] = Rq4[i][m]/NbCount[i];
			Aveq4mI[m] = Iq4[i][m]/NbCount[i];
			for(j=0;j<Nearest;j++){
				if(Nblist[i][j]<0) break;		
				Aveq4mR[m] += Rq4[Nblist[i][j]][m]/NbCount[Nblist[i][j]];
				Aveq4mI[m] += Iq4[Nblist[i][j]][m]/NbCount[Nblist[i][j]];
		}}
		Aveq6mR[m] = Rq6[i][m]/NbCount[i];
		Aveq6mI[m] = Iq6[i][m]/NbCount[i];
		for(j=0;j<Nearest;j++){
			if(Nblist[i][j]<0) break;		
			Aveq6mR[m] += Rq6[Nblist[i][j]][m]/NbCount[Nblist[i][j]];
			Aveq6mI[m] += Iq6[Nblist[i][j]][m]/NbCount[Nblist[i][j]];
	}}
	for(m=0;m<13;m++){
		if(m<9)Aveq4sq += Aveq4mR[m]*Aveq4mR[m]+Aveq4mI[m]*Aveq4mI[m];
		Aveq6sq += Aveq6mR[m]*Aveq6mR[m]+Aveq6mI[m]*Aveq6mI[m];
	}	
	for(m1=0;m1<13;m1++)for(m2=0;m2<13;m2++)for(m3=0;m3<13;m3++){
		if(m1+m2+m3==12 && m1<9 && m2<9 && m3<9){
			w4top += F4[m1][m2][m3] * (Aveq4mR[m1]*Aveq4mR[m2]*Aveq4mR[m3] -\
				Aveq4mI[m1]*Aveq4mI[m2]*Aveq4mR[m3]- Aveq4mR[m1]*Aveq4mI[m2]*Aveq4mI[m3] -\
				Aveq4mI[m1]*Aveq4mR[m2]*Aveq4mI[m3]);//Real part
//			w4topi += F4[m1][m2][m3] * (-Aveq4mR[m1]*Aveq4mR[m2]*Aveq4mI[m3] + \
				Aveq4mI[m1]*Aveq4mI[m2]*Aveq4mI[m3] + Aveq4mR[m1]*Aveq4mI[m2]*Aveq4mR[m3] +\
				Aveq4mI[m1]*Aveq4mR[m2]*Aveq4mR[m3]);//image part
		}
		if(m1+m2+m3==18){
			w6top += F6[m1][m2][m3] * (Aveq6mR[m1]*Aveq6mR[m2]*Aveq6mR[m3] - \
				Aveq6mI[m1]*Aveq6mI[m2]*Aveq6mR[m3] - Aveq6mR[m1]*Aveq6mI[m2]*Aveq6mI[m3] - \
				Aveq6mI[m1]*Aveq6mR[m2]*Aveq6mI[m3]);//Real part
//			w6topi += F6[m1][m2][m3] * (-Aveq6mR[m1]*Aveq6mR[m2]*Aveq6mI[m3] + \
				Aveq6mI[m1]*Aveq6mI[m2]*Aveq6mI[m3] + Aveq6mR[m1]*Aveq6mI[m2]*Aveq6mR[m3] + \
				Aveq6mI[m1]*Aveq6mR[m2]*Aveq6mR[m3]);//image part
	}}
	Aveq4[i] = sqrt(4.0*M_PI/9.0*Aveq4sq)/(NbCount[i]+1.0);
	Aveq6[i] = sqrt(4.0*M_PI/13.0*Aveq6sq)/(NbCount[i]+1.0);
	Avew4[i] = w4top/pow(Aveq4sq,1.5);
	Avew6[i] = w6top/pow(Aveq6sq,1.5);
}	
