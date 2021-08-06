void New_R( int, double **, double **, double *, double, double);
int Check_fail( int, int, double **, double **, double *);
double *Array1st(int );
double **Array2nd(int, int);
double ***Array3rd(int, int ,int );
double Rand();
/**********************************************************************************************************
 New_R()
 Function : 
***********************************************************************************************************/
void New_R(int N, double **ArrR_0, double **ArrR, double *L, double noise_1, double d0_1){
	int i ,k ;
	double sigma,Z;
	sigma = d0_1 * noise_1;
for(i=0;i<N;i++){
		do{	
			for(k=0;k<3;k++){
				Z = Rand();
				*(*(ArrR+i)+k) = Z * sigma + *(*(ArrR_0+i)+k);
				while( *(*(ArrR+i)+k) > L[k] ) *(*(ArrR+i)+k) -= L[k]*2.0 ;
				while( *(*(ArrR+i)+k) <-L[k] ) *(*(ArrR+i)+k) += L[k]*2.0 ;
		}}
		while(Check_fail(i,N,ArrR_0,ArrR,L));
	}
}
/**********************************************************************************************************
 New_R()
 Function : 
***********************************************************************************************************/
double Rand(){
	static int phase=1;
	static double randA, randB,Z;
	if(phase){
		randA = drand48();
		randB = drand48();
		Z = sqrt(-2.0*log(randA))*sin(2.0* M_PI* randB);
	}
	else{
		Z = sqrt(-2.0*log(randA))*cos(2.0* M_PI* randB);
	}
	phase = 1-phase;
	return Z;
}
/************************************************************************************************************
 Check_fail(I (particle number have get),Position_bin,Box_halflength)
 Function:
 Return: 1 (dis<1 redo init_R) or 0 ()
**************************************************************************************************************/
int Check_fail( int I, int N, double **ArrR_0, double **ArrR_1, double *L){
	int i, j, k;
	double diI ,dis ;
	for(i=0;i<N;i++){
		if(i==I) continue;
		diI = 0.0 ;
		if(i<I) for(k=0;k<3;k++){
			dis = *(*(ArrR_1+I)+k) - *(*(ArrR_1+i)+k); 
			while(dis > L[k]) dis -= L[k]*2.0;
			while(dis <-L[k]) dis += L[k]*2.0;
			diI += (dis*dis) ;
		}
		else for(k=0;k<3;k++){
			dis = *(*(ArrR_1+I)+k) - *(*(ArrR_0+i)+k); 
			while(dis > L[k]) dis -= L[k]*2.0;
			while(dis <-L[k]) dis += L[k]*2.0;
			diI += (dis*dis) ;
		}
		if(diI<0.3){
			return 1 ;
		}
	}
	return 0 ;
}
double *Array1st(int col){
	double *Name;
	Name = (double *)calloc(col,sizeof(double));
	return Name;
}
double **Array2nd(int col, int row){
	int i;
	double **Name;
	Name = (double **)calloc(col,sizeof(double *));
	for(i=0;i<col;i++)Name[i] = (double *)calloc(row,sizeof(double));
	return Name;
}
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

