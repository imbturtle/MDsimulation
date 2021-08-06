int ReadNLAM(FILE*);
double *ReadDataLAM(FILE*,double**);
void CalNb(int,int,int**,double**,double**,double*);
int PrintNeibor(FILE*,int,int,int**,double**,double**,double*);
void Reset(int,int**,double**);
double *Array1st(int);
double **Array2nd(int, int);
double ***Array3rd(int, int, int);
/***********************************************************************************************************
***********************************************************************************************************/
int ReadNLAM(FILE*FpR){
	char garba[100];
	int N;
	fscanf(FpR,"%[^\n] ",garba);
	fscanf(FpR,"%[^\n] ",garba);
	fscanf(FpR,"%[^\n] ",garba);
	fscanf(FpR,"%[^\n] ",garba);
	N = atoi(garba);
	rewind(FpR);
	return N;
}
/***********************************************************************************************************
***********************************************************************************************************/
double *ReadDataLAM(FILE*FpR,double**RArrR){
	char garba[100];
	int i,N;
	double RI[3],RF[3];
	static double L[3];
	fscanf(FpR,"%[^\n] ",garba);
	fscanf(FpR,"%[^\n] ",garba);
	fscanf(FpR,"%[^\n] ",garba);
	fscanf(FpR,"%[^\n] ",garba);
	N = atoi(garba);
	fscanf(FpR,"%[^\n] ",garba);
	fscanf(FpR,"%[^ ] ",garba);
	RI[0] = atof(garba);
	fscanf(FpR,"%[^\n] ",garba);
	RF[0] = atof(garba);
	fscanf(FpR,"%[^ ] ",garba);
	RI[1] = atof(garba);
	fscanf(FpR,"%[^\n] ",garba);
	RF[1] = atof(garba);
	fscanf(FpR,"%[^ ] ",garba);
	RI[2] = atof(garba);
	fscanf(FpR,"%[^\n] ",garba);
	RF[2] = atof(garba);	
	fscanf(FpR,"%[^\n] ",garba);
	for(i=0;i<3;i++) L[i] = fabs(RF[i]-RI[i])/2.;
	for(i=0;i<N;i++) fscanf(FpR,"%lf %lf %lf ",&*(*(RArrR+i)+0),&*(*(RArrR+i)+1),&*(*(RArrR+i)+2));
	return L;
}
/***********************************************************************************************************
***********************************************************************************************************/
void CalNb(int N,int i,int**Nblist,double**RArrR,double**Nbdist,double*L){
	int j,k,l;
	double Dij[3],Lij;
	extern double Rsearch;
	extern int Nsearch;
	for(j=0;j<N;j++){
		if(i==j) continue;
		for(l=0;l<3;l++){
			Dij[l] =  *(*(RArrR+j)+l) - *(*(RArrR+i)+l);
			while(Dij[l] > L[l]) Dij[l] -= L[l]*2.0;
			while(Dij[l] <-L[l]) Dij[l] += L[l]*2.0;
		}
		Lij = sqrt(Dij[0] * Dij[0] + Dij[1] * Dij[1] + Dij[2] * Dij[2]);
		if(Lij > Rsearch) continue;
		for(k=0;k<Nsearch;k++)if(Lij < Nbdist[i][k]){
			for(l=Nsearch;l>k;l--){
				Nbdist[i][l] = Nbdist[i][l-1];
				Nblist[i][l] = Nblist[i][l-1];
			}
			Nbdist[i][k] = Lij;
			Nblist[i][k] = j;
			break;
}}}
/***********************************************************************************************************
***********************************************************************************************************/
int PrintNeibor(FILE*FpW,int N,int i,int**Nblist,double**RArrR,double**Nbdist,double*L){
	int j,l,Nbcount=0;
	double Dij[3];
	extern int Nsearch,Ncut;
	extern double Rsearch,Rcut;
	for(j=0;j<Nsearch;j++){
		if(Nblist[i][j]>=0){
			for(l=0;l<3;l++){
				Dij[l] =  RArrR[Nblist[i][j]][l] - RArrR[i][l];
				while(Dij[l] > L[l]) Dij[l] -= L[l]*2.0;
				while(Dij[l] <-L[l]) Dij[l] += L[l]*2.0;
			}
			if(Nbdist[i][j]/Nbdist[i][0]<=Rcut && j<Ncut) fprintf(FpW,"%lf %lf %lf\n",Dij[0]/Nbdist[i][0],Dij[1]/Nbdist[i][0],Dij[2]/Nbdist[i][0]);
			if(Nbdist[i][j]/Nbdist[i][0]<=Rcut) Nbcount++;
			if(Nbdist[i][j]/Nbdist[i][0]>Rcut && j<Ncut) fprintf(FpW,"%lf %lf %lf\n",0.0,0.0,0.0);
		}
		if(Nblist[i][j]<0 && j<Ncut) fprintf(FpW,"%lf %lf %lf\n",0.0,0.0,0.0);
	}
	return Nbcount;
}
/***********************************************************************************************************
***********************************************************************************************************/
void Reset(int N,int **Nblist,double **Nbdist){
	int i,j;
	extern int Nsearch;
	for(i=0;i<N;i++)for(j=0;j<Nsearch;j++){
		Nbdist[i][j] = 99.0;
		Nblist[i][j] = -1;
}}
/***********************************************************************************************************
***********************************************************************************************************/
double *Array1st(int col){
	double *Name;
	Name = (double *)calloc(col,sizeof(double));
	return Name;
}
/***********************************************************************************************************
***********************************************************************************************************/
double **Array2nd(int col, int row){
	int i;
	double **Name;
	Name = (double **)calloc(col,sizeof(double *));
	for(i=0;i<col;i++)Name[i] = (double *)calloc(row,sizeof(double));
	return Name;
}
/***********************************************************************************************************
***********************************************************************************************************/
double ***Array3rd(int col, int row,int hight){
	int i,j;
	double ***Name;
	Name = (double ***)calloc(col,sizeof(double **));
	for(i=0;i<col;i++){
		Name[i] = (double **)calloc(row,sizeof(double *));
		for(j=0;j<row;j++) Name[i][j] = (double *)calloc(hight,sizeof(double));
	}
	return Name;
}
/***********************************************************************************************************
***********************************************************************************************************/
