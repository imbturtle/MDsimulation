#include<stdlib.h>
#include<math.h>
#include"Xcate_Ave.h"
/***********************************************************************************************************
***********************************************************************************************************/
void CrystalAvevec(int N,int m,double*halfL,int**Nblist,double**ArrR,double**LocVecq4R,double**LocVecq4I,double**LocVecq6R,double**LocVecq6I\
	,double*eigenvalue,double**matric,double**Xcate,double*NbCount){
	int i,j,lt1;
	double a,b,c,cosa,cosb,cosc,sina,sinb,sinc,normal,testa,testc;	
	double *iv,*niv,*jv,*njv,*kv,*nkv,*n,*Rmatrix;
	double **NArrR,XS=0.0;
	
	Rmatrix = Array1st(9);
	iv = Array1st(3);
	niv = Array1st(3);
	jv = Array1st(3);
	njv = Array1st(3);
	kv = Array1st(3);
	nkv = Array1st(3);
	n = Array1st(3);
	NArrR = Array2nd(N,3);
	for(i=1;i<4;i++)for(j=1;j<4;j++)if(eigenvalue[i]>eigenvalue[j]){
		a=eigenvalue[j];
		eigenvalue[j]=eigenvalue[i];
		eigenvalue[i]=a;
		a=matric[1][j];
		matric[1][j]=matric[1][i];
		matric[1][i]=a;
		a=matric[2][j];
		matric[2][j]=matric[2][i];
		matric[2][i]=a;
		a=matric[3][j];
		matric[3][j]=matric[3][i];
		matric[3][i]=a;
	}
	for(i=0;i<3;i++) iv[i]=jv[i]=kv[i]=0.0;
	iv[0]=jv[1]=kv[2]=1.0;
	niv[0]=matric[1][1];
	niv[1]=matric[2][1];
	niv[2]=matric[3][1];
	njv[0]=matric[1][2];
	njv[1]=matric[2][2];
	njv[2]=matric[3][2];
	nkv[0]=matric[2][1]*matric[3][2]-matric[3][1]*matric[2][2];
	nkv[1]=matric[3][1]*matric[1][2]-matric[1][1]*matric[3][2];
	nkv[2]=matric[1][1]*matric[2][2]-matric[2][1]*matric[1][2];
	n[0]=kv[1]*nkv[2]-kv[2]*nkv[1];
	n[1]=kv[2]*nkv[0]-kv[0]*nkv[2];
	n[2]=kv[0]*nkv[1]-kv[1]*nkv[0];
	normal=pow((n[0]*n[0]+n[1]*n[1]+n[2]*n[2]),0.5);
	n[0]=n[0]/normal;
	n[1]=n[1]/normal;
	n[2]=n[2]/normal;
	if (normal==0.0) n[0]=n[1]=n[2]=0.0;
	cosb=kv[0]*nkv[0]+kv[1]*nkv[1]+kv[2]*nkv[2];
	cosa=n[0]*iv[0]+n[1]*iv[1]+n[2]*iv[2];
	cosc=n[0]*niv[0]+n[1]*niv[1]+n[2]*niv[2];
	testa=0.0;
	for(i=0;i<3;i++) testa+=n[i]*jv[i];
	a=acos(cosa)*180.0/M_PI;
	if(testa<0.0) a=a+180.0;
    b=acos(cosb)*180.0/M_PI;
    c=acos(cosc)*180.0/M_PI;
	testc=0.0;
	for(i=0;i<3;i++) testc+=n[i]*njv[i];
	if(testc>0.0) c=c+180.0;
	sina=sin(a*M_PI/180.0);
	sinb=sin(b*M_PI/180.0);
	sinc=sin(c*M_PI/180.0);
	Rmatrix[0]=cosc*cosa-cosb*sina*sinc;
	Rmatrix[3]=-sinc*cosa-cosb*sina*cosc;
	Rmatrix[6]=sinb*sina;
	Rmatrix[1]=cosc*sina+cosb*cosa*sinc;
	Rmatrix[4]=-sinc*sina+cosa*cosb*cosc;
	Rmatrix[7]=-sinb*cosa;
	Rmatrix[2]=sinb*sinc;
	Rmatrix[5]=sinb*cosc;
	Rmatrix[8]=cosb;
	for(i=0;i<N;i++){
		NArrR[i][0]=Rmatrix[0]*ArrR[i][0]+Rmatrix[1]*ArrR[i][1]+Rmatrix[2]*ArrR[i][2];
		NArrR[i][1]=Rmatrix[3]*ArrR[i][0]+Rmatrix[4]*ArrR[i][1]+Rmatrix[5]*ArrR[i][2];
		NArrR[i][2]=Rmatrix[6]*ArrR[i][0]+Rmatrix[7]*ArrR[i][1]+Rmatrix[8]*ArrR[i][2];
	}
	Cal_AveVecQ(N,m,Nblist,NArrR,LocVecq4R,LocVecq4I,LocVecq6R,LocVecq6I,halfL,Rmatrix,NbCount);
	LSS(m,LocVecq4R,LocVecq4I,LocVecq6R,LocVecq6I,Xcate);
	free(iv);
	free(niv);
	free(jv);
	free(njv);
	free(kv);
	free(nkv);
	free(n);
	free(Rmatrix);
	for(j=0;j<N;j++){
		free(NArrR[j]);
	}
	free(NArrR);
}
