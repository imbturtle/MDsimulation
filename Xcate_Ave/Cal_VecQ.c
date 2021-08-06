#include<stdlib.h>
#include<math.h>
#include"Xcate_Ave.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Cal_VecQ(int N,int i,int**Nblist,double**ArrR,double**LocVecq4R,double**LocVecq4I,double**LocVecq6R,double**LocVecq6I,double*L,double*Rmatrix){
	int j,l,m;
	double Lij,Dij[3],costheta,sintheta,phi;
	double **Rq4,**Iq4,**Rq6,**Iq6,*Locq4sq,*Locq6sq;
	extern int Nearest;
	Rq4 = Array2nd(N,9);
	Iq4 = Array2nd(N,9);
	Rq6 = Array2nd(N,13);
	Iq6 = Array2nd(N,13);
	Locq4sq = Array1st(N);
	Locq6sq = Array1st(N);
	for(j=0;j<Nearest;j++){
		if(Nblist[i][j]<0) break;
		for(l=0;l<3;l++)Dij[l] =  ArrR[Nblist[i][j]][l] - ArrR[i][l];
		while(Dij[0]*Rmatrix[0]+Dij[1]*Rmatrix[3]+Dij[2]*Rmatrix[6] > L[0]){
			Dij[0] -= (Rmatrix[0]*2.0*L[0]);
			Dij[1] -= (Rmatrix[3]*2.0*L[0]);
			Dij[2] -= (Rmatrix[6]*2.0*L[0]);
		}
		while(Dij[0]*Rmatrix[0]+Dij[1]*Rmatrix[3]+Dij[2]*Rmatrix[6] < -L[0]){
			Dij[0] += (Rmatrix[0]*2.0*L[0]);
			Dij[1] += (Rmatrix[3]*2.0*L[0]);
			Dij[2] += (Rmatrix[6]*2.0*L[0]);
		}
		while(Dij[0]*Rmatrix[1]+Dij[1]*Rmatrix[4]+Dij[2]*Rmatrix[7] > L[1]){
			Dij[0] -= (Rmatrix[1]*2.0*L[1]);
			Dij[1] -= (Rmatrix[4]*2.0*L[1]);
			Dij[2] -= (Rmatrix[7]*2.0*L[1]);
		}
		while(Dij[0]*Rmatrix[1]+Dij[1]*Rmatrix[4]+Dij[2]*Rmatrix[7] < -L[1]){
			Dij[0] += (Rmatrix[1]*2.0*L[1]);
			Dij[1] += (Rmatrix[4]*2.0*L[1]);
			Dij[2] += (Rmatrix[7]*2.0*L[1]);
		}
		while(Dij[0]*Rmatrix[2]+Dij[1]*Rmatrix[5]+Dij[2]*Rmatrix[8] > L[2]){
			Dij[0] -= (Rmatrix[2]*2.0*L[2]);
			Dij[1] -= (Rmatrix[5]*2.0*L[2]);
			Dij[2] -= (Rmatrix[8]*2.0*L[2]);
		}
		while(Dij[0]*Rmatrix[2]+Dij[1]*Rmatrix[5]+Dij[2]*Rmatrix[8] < -L[2]){
			Dij[0] += (Rmatrix[2]*2.0*L[2]);
			Dij[1] += (Rmatrix[5]*2.0*L[2]);
			Dij[2] += (Rmatrix[8]*2.0*L[2]);
		}
		Lij = sqrt(Dij[0] * Dij[0] + Dij[1] * Dij[1] + Dij[2] * Dij[2]);
		costheta = Dij[2]/Lij;
		sintheta = pow(1.0-costheta*costheta,0.5);
		phi = atan(Dij[1]/Dij[0]);
		if(Dij[0]==0.0 && Dij[1]==0.0) phi = 0.0;
		if(Dij[0]==0.0 && Dij[1]>0.0)  phi = M_PI/2.0;
		if(Dij[0]==0.0 && Dij[1]<0.0)  phi =-M_PI/2.0;
		if(phi < 0.0) phi += 2.0*M_PI;
		if(Dij[0] < 0.0) phi += M_PI;
		Rq4[i][0] += (3.0/16.0)*sqrt(35.0/2.0/M_PI)*pow(sintheta,4.0)*cos(4.0*phi);
		Iq4[i][0] += (3.0/16.0)*sqrt(35.0/2.0/M_PI)*pow(sintheta,4.0)*-sin(4.0*phi);
		Rq4[i][1] += (3.0/8.0)*sqrt(35.0/M_PI)*pow(sintheta,3.0)*costheta*cos(3.0*phi);
		Iq4[i][1] += (3.0/8.0)*sqrt(35.0/M_PI)*pow(sintheta,3.0)*costheta*-sin(3.0*phi);
		Rq4[i][2] += (3.0/8.0)*sqrt(5.0/2.0/M_PI)*pow(sintheta,2.0)*(7.0*pow(costheta,2.0)-1)*cos(2.0*phi);
		Iq4[i][2] += (3.0/8.0)*sqrt(5.0/2.0/M_PI)*pow(sintheta,2.0)*(7.0*pow(costheta,2.0)-1)*-sin(2.0*phi);
		Rq4[i][3] += (3.0/8.0)*sqrt(5.0/M_PI)*sintheta*(7.0*pow(costheta,3.0)-3.0*costheta)*cos(phi);
		Iq4[i][3] += (3.0/8.0)*sqrt(5.0/M_PI)*sintheta*(7.0*pow(costheta,3.0)-3.0*costheta)*-sin(phi);
		Rq4[i][4] += (3.0/16.0)*sqrt(1.0/M_PI)*(35.0*pow(costheta,4.0)-30.0*pow(costheta,2.0)+3.0);
		Rq4[i][5] += (-3.0/8.0)*sqrt(5.0/M_PI)*sintheta*(7.0*pow(costheta,3.0)-3.0*costheta)*cos(phi);
		Iq4[i][5] += (-3.0/8.0)*sqrt(5.0/M_PI)*sintheta*(7.0*pow(costheta,3.0)-3.0*costheta)*sin(phi);
		Rq4[i][6] += (3.0/8.0)*sqrt(5.0/2.0/M_PI)*pow(sintheta,2.0)*(7.0*pow(costheta,2.0)-1)*cos(2.0*phi);
		Iq4[i][6] += (3.0/8.0)*sqrt(5.0/2.0/M_PI)*pow(sintheta,2.0)*(7.0*pow(costheta,2.0)-1)*sin(2.0*phi);
		Rq4[i][7] += (-3.0/8.0)*sqrt(35.0/M_PI)*pow(sintheta,3.0)*costheta*cos(3.0*phi);
		Iq4[i][7] += (-3.0/8.0)*sqrt(35.0/M_PI)*pow(sintheta,3.0)*costheta*sin(3.0*phi);
		Rq4[i][8] += (3.0/16.0)*sqrt(35.0/2.0/M_PI)*pow(sintheta,4.0)*cos(4.0*phi);
		Iq4[i][8] += (3.0/16.0)*sqrt(35.0/2.0/M_PI)*pow(sintheta,4.0)*sin(4.0*phi);
		Rq6[i][0] += (1.0/64.0)*sqrt(3003.0/M_PI)*pow(sintheta,6.0)*cos(6.0*phi);
		Iq6[i][0] += (1.0/64.0)*sqrt(3003.0/M_PI)*pow(sintheta,6.0)*-sin(6.0*phi);
		Rq6[i][1] += (3.0/32.0)*sqrt(1001.0/M_PI)*pow(sintheta,5.0)*costheta*cos(5.0*phi);
		Iq6[i][1] += (3.0/32.0)*sqrt(1001.0/M_PI)*pow(sintheta,5.0)*costheta*-sin(5.0*phi);
		Rq6[i][2] += (3.0/32.0)*sqrt(91.0/2.0/M_PI)*pow(sintheta,4.0)*(11.0*pow(costheta,2.0)-1.0)*cos(4.0*phi);
		Iq6[i][2] += (3.0/32.0)*sqrt(91.0/2.0/M_PI)*pow(sintheta,4.0)*(11.0*pow(costheta,2.0)-1.0)*-sin(4.0*phi);
		Rq6[i][3] += (1.0/32.0)*sqrt(1365.0/M_PI)*pow(sintheta,3.0)*(11.0*pow(costheta,3.0)-3.0*costheta)*cos(3.0*phi);
		Iq6[i][3] += (1.0/32.0)*sqrt(1365.0/M_PI)*pow(sintheta,3.0)*(11.0*pow(costheta,3.0)-3.0*costheta)*-sin(3.0*phi);
		Rq6[i][4] += (1.0/64.0)*sqrt(1365.0/M_PI)*pow(sintheta,2.0)*\
			(33.0*pow(costheta,4.0)-18.0*pow(costheta,2.0)+1.0)*cos(2.0*phi);
		Iq6[i][4] += (1.0/64.0)*sqrt(1365.0/M_PI)*pow(sintheta,2.0)*\
			(33.0*pow(costheta,4.0)-18.0*pow(costheta,2.0)+1.0)*-sin(2.0*phi);
		Rq6[i][5] += (1.0/16.0)*sqrt(273.0/2.0/M_PI)*sintheta*\
			(33.0*pow(costheta,5.0)-30.0*pow(costheta,3.0)+5.0*costheta)*cos(phi);
		Iq6[i][5] += (1.0/16.0)*sqrt(273.0/2.0/M_PI)*sintheta*\
			(33.0*pow(costheta,5.0)-30.0*pow(costheta,3.0)+5.0*costheta)*-sin(phi);
		Rq6[i][6] += (1.0/32.0)*sqrt(13.0/M_PI)*(231.0*pow(costheta,6.0)-315.0*pow(costheta,4.0)+\
			105.0*pow(costheta,2.0)-5);
		Rq6[i][7] += (-1.0/16.0)*sqrt(273.0/2.0/M_PI)*sintheta*\
			(33.0*pow(costheta,5.0)-30.0*pow(costheta,3.0)+5.0*costheta)*cos(phi);
		Iq6[i][7] += (-1.0/16.0)*sqrt(273.0/2.0/M_PI)*sintheta*\
			(33.0*pow(costheta,5.0)-30.0*pow(costheta,3.0)+5.0*costheta)*sin(phi);
		Rq6[i][8] += (1.0/64.0)*sqrt(1365.0/M_PI)*pow(sintheta,2.0)*\
			(33.0*pow(costheta,4.0)-18.0*pow(costheta,2.0)+1.0)*cos(2.0*phi);
		Iq6[i][8] += (1.0/64.0)*sqrt(1365.0/M_PI)*pow(sintheta,2.0)*\
			(33.0*pow(costheta,4.0)-18.0*pow(costheta,2.0)+1.0)*sin(2.0*phi);
		Rq6[i][9] += (-1.0/32.0)*sqrt(1365.0/M_PI)*pow(sintheta,3.0)*(11.0*pow(costheta,3.0)-3.0*costheta)*cos(3.0*phi);
		Iq6[i][9] += (-1.0/32.0)*sqrt(1365.0/M_PI)*pow(sintheta,3.0)*(11.0*pow(costheta,3.0)-3.0*costheta)*sin(3.0*phi);
		Rq6[i][10] += (3.0/32.0)*sqrt(91.0/2.0/M_PI)*pow(sintheta,4.0)*(11.0*pow(costheta,2.0)-1.0)*cos(4.0*phi);
		Iq6[i][10] += (3.0/32.0)*sqrt(91.0/2.0/M_PI)*pow(sintheta,4.0)*(11.0*pow(costheta,2.0)-1.0)*sin(4.0*phi);
		Rq6[i][11] += (-3.0/32.0)*sqrt(1001.0/M_PI)*pow(sintheta,5.0)*costheta*cos(5.0*phi);
		Iq6[i][11] += (-3.0/32.0)*sqrt(1001.0/M_PI)*pow(sintheta,5.0)*costheta*sin(5.0*phi);
		Rq6[i][12] += (1.0/64.0)*sqrt(3003.0/M_PI)*pow(sintheta,6.0)*cos(6.0*phi);
		Iq6[i][12] += (1.0/64.0)*sqrt(3003.0/M_PI)*pow(sintheta,6.0)*sin(6.0*phi); 						
	}
	for(m=0;m<9;m++)Locq4sq[i] += Rq4[i][m]*Rq4[i][m]+Iq4[i][m]*Iq4[i][m];
	for(m=0;m<13;m++)Locq6sq[i] += Rq6[i][m]*Rq6[i][m]+Iq6[i][m]*Iq6[i][m];
	for(m=0;m<13;m++){
		if(m<9){
			LocVecq4R[i][m] = Rq4[i][m]/sqrt(Locq4sq[i]);
			LocVecq4I[i][m] = Iq4[i][m]/sqrt(Locq4sq[i]);
		}
		LocVecq6R[i][m] = Rq6[i][m]/sqrt(Locq6sq[i]);
		LocVecq6I[i][m] = Iq6[i][m]/sqrt(Locq6sq[i]);
	}
	for(j=0;j<N;j++){
		free(Rq4[j]);
		free(Iq4[j]);
		free(Rq6[j]);
		free(Iq6[j]);
	}
	free(Rq4);
	free(Iq4);
	free(Rq6);
	free(Iq6);
	free(Locq4sq);
	free(Locq6sq);
}
