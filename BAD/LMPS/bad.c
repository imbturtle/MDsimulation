int ReadN(FILE*);
double *ReadData(FILE*,double**);
void CalBAD(int,int,double*,double**,double*,double*);
double *Array1st(int);
double **Array2nd(int,int);
double ***Array3rd(int,int,int);
#define pi M_PI
/***********************************************************************************************************
***********************************************************************************************************/
int ReadN(FILE*FpR){
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
double *ReadData(FILE*FpR,double**RArrR){
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
void CalBAD(int N,int i,double*L,double**RArrR,double*Angle,double*CN){
	int j,k,l,Count=0;
	double Dij[3],Dik[3],Lij,Lik,Rad=180.0/M_PI,TotalCount=0.,NAngle;
	extern double Rmin,Dtheta;
	for(j=i+1;j<N;j++){//j top
		if(i==j) continue;
		for(l=0;l<3;l++){
			Dij[l] =  RArrR[j][l] - RArrR[i][l];
			while(Dij[l] > L[l]) Dij[l] -= L[l]*2.0;
			while(Dij[l] <-L[l]) Dij[l] += L[l]*2.0;
		}
		Lij = sqrt(Dij[0] * Dij[0] + Dij[1] * Dij[1] + Dij[2] * Dij[2]);
		if(Lij > Rmin) continue;
		Count++;
		for(k=j+1;k<N;k++){//k top
			if(i==k) continue;
			for(l=0;l<3;l++){
				Dik[l] =  RArrR[k][l] - RArrR[i][l];
				while(Dik[l] > L[l]) Dik[l] -= L[l]*2.0;
				while(Dik[l] <-L[l]) Dik[l] += L[l]*2.0;
			}
			Lik = sqrt(Dik[0] * Dik[0] + Dik[1] * Dik[1] + Dik[2] * Dik[2]);  
			if(Lik > Rmin) continue; 
			NAngle=((acos((Dij[0] * Dik[0] + Dij[1] * Dik[1] + Dij[2] * Dik[2])/(Lij * Lik)) * Rad)/Dtheta);
			if(((Dij[0] * Dik[0] + Dij[1] * Dik[1] + Dij[2] * Dik[2])/(Lij * Lik))<-1.0) NAngle=179.9;
			Angle[(int)NAngle]++;
	}}
	CN[Count]++;
}
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
