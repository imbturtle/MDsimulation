#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"Xcate_Ave.h"
/***********************************************************************************************************
***********************************************************************************************************/
void LSS(int i,double**LocVecq4R,double**LocVecq4I,double**LocVecq6R,double**LocVecq6I,double**Xcate){
	int j,l,m,**NbBCClist,**NbFCClist,**NbHCPlist,NBCC=15,NFCC=13,NHCP=13;
	double **ArrRBCC,**ArrRFCC,**ArrRHCP,**NArrRBCC,**NArrRFCC,**NArrRHCP,*Rmatrix\
	,**BCCLocVecq4R,**BCCLocVecq4I,**BCCLocVecq6R,**BCCLocVecq6I\
	,**FCCLocVecq4R,**FCCLocVecq4I,**FCCLocVecq6R,**FCCLocVecq6I\
	,**HCPLocVecq4R,**HCPLocVecq4I,**HCPLocVecq6R,**HCPLocVecq6I;
	double a,b,c,cosa,cosb,cosc,sina,sinb,sinc,X4,X6,L[3]={99.0,99.0,99.0},Rot_interval=5.0;
	extern int Nearest;
	FILE *FpRcry;
	Rmatrix = Array1st(9);
	FpRcry = fopen("./Crystal/BCC.txt","r");
	ArrRBCC = Array2nd(NBCC,3);
	NArrRBCC = Array2nd(NBCC,3);
	NbBCClist = (int **)calloc(NBCC,sizeof(int *));
	for(j=0;j<NBCC;j++) NbBCClist[j] = (int *)calloc(Nearest,sizeof(int));
	BCCLocVecq4R = Array2nd(NBCC,9);
	BCCLocVecq4I = Array2nd(NBCC,9);
	BCCLocVecq6R = Array2nd(NBCC,13);
	BCCLocVecq6I = Array2nd(NBCC,13);
	for(j=0;j<NBCC;j++) for(l=0;l<Nearest;l++) NbBCClist[j][l]=-1;
	for(j=0;j<NBCC;j++)fscanf(FpRcry,"%lf %lf %lf",&*(*(ArrRBCC+j)+0),&*(*(ArrRBCC+j)+1),&*(*(ArrRBCC+j)+2));
	fclose(FpRcry);
	FpRcry = fopen("./Crystal/FCC.txt","r");
	ArrRFCC = Array2nd(NFCC,3);
	NArrRFCC = Array2nd(NFCC,3);
	NbFCClist = (int **)calloc(NFCC,sizeof(int *));
	for(j=0;j<NFCC;j++) NbFCClist[j] = (int *)calloc(Nearest,sizeof(int));
	FCCLocVecq4R = Array2nd(NFCC,9);
	FCCLocVecq4I = Array2nd(NFCC,9);
	FCCLocVecq6R = Array2nd(NFCC,13);
	FCCLocVecq6I = Array2nd(NFCC,13);
	for(j=0;j<NFCC;j++) for(l=0;l<Nearest;l++) NbFCClist[j][l]=-1;
	for(j=0;j<NFCC;j++)fscanf(FpRcry,"%lf %lf %lf",&*(*(ArrRFCC+j)+0),&*(*(ArrRFCC+j)+1),&*(*(ArrRFCC+j)+2));
	fclose(FpRcry);
	FpRcry = fopen("./Crystal/HCP.txt","r");
	ArrRHCP = Array2nd(13,3);
	NArrRHCP = Array2nd(NHCP,3);
	NbHCPlist = (int **)calloc(NHCP,sizeof(int *));
	for(j=0;j<NHCP;j++) NbHCPlist[j] = (int *)calloc(Nearest,sizeof(int));
	HCPLocVecq4R = Array2nd(NHCP,9);
	HCPLocVecq4I = Array2nd(NHCP,9);
	HCPLocVecq6R = Array2nd(NHCP,13);
	HCPLocVecq6I = Array2nd(NHCP,13);
	for(j=0;j<NHCP;j++) for(l=0;l<Nearest;l++) NbHCPlist[j][l]=-1;
	for(j=0;j<NHCP;j++)fscanf(FpRcry,"%lf %lf %lf",&*(*(ArrRHCP+j)+0),&*(*(ArrRHCP+j)+1),&*(*(ArrRHCP+j)+2));
	fclose(FpRcry);
	for(a=0.0;a<360.0;a+=Rot_interval){
		for(b=0.0;b<180.0;b+=Rot_interval){
			for(c=0.0;c<360.0;c+=Rot_interval){
				cosa=cos(a*M_PI/180.0);
				cosb=cos(b*M_PI/180.0);
				cosc=cos(c*M_PI/180.0);
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
				for(j=0;j<NBCC;j++) for(l=0;l<3;l++) NArrRBCC[j][l]=Rmatrix[0+l*3]*ArrRBCC[j][0]+Rmatrix[1+l*3]*ArrRBCC[j][1]+Rmatrix[2+l*3]*ArrRBCC[j][2];
				Cal_CryNb(NBCC,0,NbBCClist,NArrRBCC,L,9.0);
				Cal_VecQ(NBCC,0,NbBCClist,NArrRBCC,BCCLocVecq4R,BCCLocVecq4I,BCCLocVecq6R,BCCLocVecq6I,L,Rmatrix);
				X4=0.0;
				X6=0.0;
				for(m=0;m<9;m++) X4 += LocVecq4R[i][m]*BCCLocVecq4R[0][m] + LocVecq4I[i][m]*BCCLocVecq4I[0][m];
				for(m=0;m<13;m++) X6 += LocVecq6R[i][m]*BCCLocVecq6R[0][m] + LocVecq6I[i][m]*BCCLocVecq6I[0][m];
				if((X4>0.0)&&(X6>0.0)&&(X4*X6>Xcate[i][0])) Xcate[i][0]=X4*X6;
				for(j=0;j<NFCC;j++) for(l=0;l<3;l++) NArrRFCC[j][l]=Rmatrix[0+l*3]*ArrRFCC[j][0]+Rmatrix[1+l*3]*ArrRFCC[j][1]+Rmatrix[2+l*3]*ArrRFCC[j][2];
				Cal_CryNb(NFCC,0,NbFCClist,NArrRFCC,L,9.0);
				Cal_VecQ(NFCC,0,NbFCClist,NArrRFCC,FCCLocVecq4R,FCCLocVecq4I,FCCLocVecq6R,FCCLocVecq6I,L,Rmatrix);
				X4=0.0;
				X6=0.0;
				for(m=0;m<9;m++) X4 += LocVecq4R[i][m]*FCCLocVecq4R[0][m] + LocVecq4I[i][m]*FCCLocVecq4I[0][m];
				for(m=0;m<13;m++) X6 += LocVecq6R[i][m]*FCCLocVecq6R[0][m] + LocVecq6I[i][m]*FCCLocVecq6I[0][m];
				if((X4>0.0)&&(X6>0.0)&&(X4*X6>Xcate[i][1])) Xcate[i][1]=X4*X6;
				for(j=0;j<NHCP;j++) for(l=0;l<3;l++) NArrRHCP[j][l]=Rmatrix[0+l*3]*ArrRHCP[j][0]+Rmatrix[1+l*3]*ArrRHCP[j][1]+Rmatrix[2+l*3]*ArrRHCP[j][2];
				Cal_CryNb(NHCP,0,NbHCPlist,NArrRHCP,L,9.0);
				Cal_VecQ(NHCP,0,NbHCPlist,NArrRHCP,HCPLocVecq4R,HCPLocVecq4I,HCPLocVecq6R,HCPLocVecq6I,L,Rmatrix);
				X4=0.0;
				X6=0.0;
				for(m=0;m<9;m++) X4 += LocVecq4R[i][m]*HCPLocVecq4R[0][m] + LocVecq4I[i][m]*HCPLocVecq4I[0][m];
				for(m=0;m<13;m++) X6 += LocVecq6R[i][m]*HCPLocVecq6R[0][m] + LocVecq6I[i][m]*HCPLocVecq6I[0][m];
				if((X4>0.0)&&(X6>0.0)&&(X4*X6>Xcate[i][2])) Xcate[i][2]=X4*X6;
	}}}
	free(Rmatrix);
	for(j=0;j<15;j++){
		free(NbBCClist[j]);
		free(ArrRBCC[j]);
		free(NArrRBCC[j]);
		free(BCCLocVecq4R[j]);
		free(BCCLocVecq4I[j]);
		free(BCCLocVecq6I[j]);
		free(BCCLocVecq6R[j]);
	}
	free(NbBCClist);
	free(ArrRBCC);
	free(NArrRBCC);
	free(BCCLocVecq4R);
	free(BCCLocVecq4I);
	free(BCCLocVecq6I);
	free(BCCLocVecq6R);
	
	for(j=0;j<NFCC;j++){
		free(NbFCClist[j]);
		free(ArrRFCC[j]);
		free(NArrRFCC[j]);
		free(FCCLocVecq4R[j]);
		free(FCCLocVecq4I[j]);
		free(FCCLocVecq6I[j]);
		free(FCCLocVecq6R[j]);
	}
	free(NbFCClist);
	free(ArrRFCC);
	free(NArrRFCC);
	free(FCCLocVecq4R);
	free(FCCLocVecq4I);
	free(FCCLocVecq6I);
	free(FCCLocVecq6R);
	
	for(j=0;j<NHCP;j++){
		free(NbHCPlist[j]);
		free(ArrRHCP[j]);
		free(NArrRHCP[j]);
		free(HCPLocVecq4R[j]);
		free(HCPLocVecq4I[j]);
		free(HCPLocVecq6I[j]);
		free(HCPLocVecq6R[j]);
	}
	free(NbHCPlist);
	free(ArrRHCP);
	free(NArrRHCP);
	free(HCPLocVecq4R);
	free(HCPLocVecq4I);
	free(HCPLocVecq6I);
	free(HCPLocVecq6R);
}

