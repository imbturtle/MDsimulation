#include<stdio.h>
#include<math.h>
#include"MD.h"
/*****************************************************************************************************************
 halfleng(Nparticles,density)
 Function:Decide halfleng
******************************************************************************************************************/
double BoxLeng(int N,double d){
	double V ,halfL ;
	V=N/d ;
	halfL=pow(V,1.0/3.0)/2.0;
	return halfL;
}
