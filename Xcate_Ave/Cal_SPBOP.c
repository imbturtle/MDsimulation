#include<math.h>
#include"Xcate_Ave.h"
/***********************************************************************************************************
***********************************************************************************************************/
double Cal_SPBOP(int N,int i,int*Connect,int**Nblist,double**LocVecq4R,double**LocVecq4I,double**LocVecq6R,double**LocVecq6I){
	int j,m,Count=0;
	double Si;
	extern int Nearest;
	extern double qc;
	Connect[i]=0;
	for(j=0;j<Nearest;j++){
		if(Nblist[i][j]<0) break;
		Si=0.;
		Count++;
		for(m=0;m<13;m++) Si += LocVecq6R[i][m]*LocVecq6R[Nblist[i][j]][m] + LocVecq6I[i][m]*LocVecq6I[Nblist[i][j]][m];
		if(Si>qc) Connect[i]++;
	}
	Si /= Count;
	return Si;
}
