#include"Xcate_Ave.h"
/***********************************************************************************************************
***********************************************************************************************************/
double SumPi(double x){
	int i,j;
	if(x == 0.0){
		return 1.0;
	}
	j = (int)x;
	x = 1.0;
	for(i=1;i<=j;i++) x *= i;
	return x;
}
