double Cal_Nb(int,int,int**,double**,double**,double*);
void Cal_CryNb(int,int,int**,double**,double*,double);
void Cal_SHP(int,int,int**,double**,double**,double**,double**,double**,double*,double*,double**,double**,double**,double**,double*);
void Cal_LocQW(int,int,int**,double*,double**,double**,double**,double**,double***,double***\
	,double*,double*,double*,double*,double*,double*);

void Cal_VecQ(int,int,int**,double**,double**,double**,double**,double**,double*,double*);
void Cal_AveVecQ(int N,int i,int**,double**,double**,double**,double**,double**,double*,double*,double*);
double Cal_SPBOP(int,int,int*,int**,double**,double**,double**,double**);
void Cal_AveQW(int,int,int**,double*,double**,double**,double**,double**,double***,double***,double*,double*,double*,double*);
void Reset(int,int**,double**,double**,double**,double**,double**,double*,double*,double**);
void Winger3j(double***,double***);
void SumYml(double**,double**,double**,double**,double*,double,int);
double SumPi(double);
double *Array1st(int);
double **Array2nd(int,int);
double ***Array3rd(int,int,int);
void Crystalvec(int,int,double*,int**,double**,double**,double**,double**,double**,double*,double**,double**);
void LSS(int,double**,double**,double**,double**,double**);
void mat(int,int**,double*,double**,double**);
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void tqli(double *, double *, long , double **);
double pythag(double, double);
void tred2(double **, long, double *, double *);
