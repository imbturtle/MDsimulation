#include<math.h>
#include"Xcate.h"
/***********************************************************************************************************
***********************************************************************************************************/
void Winger3j(double ***F4, double ***F6){
	double diagonal, summation, m1, m2, m3, l1, l2, l3, t;
	l1=l2=l3=4.0;
	for(m1=-l1;m1<l1+1.0;m1++){
		for(m2=-l2;m2<l2+1.0;m2++){
			for(m3=-l3;m3<l3+1.0;m3++){
				if(m1+m2+m3 == 0.0){
					summation = 0.0;
					for(t=0.0;t<=20.0;t+=1.0){
						if(l3-l2+t+m1<0.0 || l3-l1+t-m2<0.0 || l1+l2-l3-t<0.0 || l1-t-m1<0.0 || l2-t+m2 <0.0) continue;
						summation += pow(-1.0,t)/(SumPi(t)*SumPi(l3-l2+t+m1)*SumPi(l3-l1+t-m2)*SumPi(l1+l2-l3-t)*\
						SumPi(l1-t-m1)*SumPi(l2-t+m2));
					}
					diagonal = pow(-1.0,(double)(l1-l2-m3)) *\
					sqrt(SumPi(l1+l2-l3)*SumPi(l1-l2+l3)*SumPi(-l1+l2+l3)/SumPi(l1+l2+l3+1.0)) *\
					sqrt(SumPi(l1+m1)*SumPi(l1-m1)*SumPi(l2+m2)*SumPi(l2-m2)*SumPi(l3+m3)*SumPi(l3-m3)) *\
					summation;
					F4[(int)(m1+l1)][(int)(m2+l2)][(int)(m3+l3)] = diagonal;
	}}}}
	l1=l2=l3=6.0;
	for(m1=-l1;m1<l1+1.0;m1++){
		for(m2=-l2;m2<l2+1.0;m2++){
			for(m3=-l3;m3<l3+1.0;m3++){
				if(m1+m2+m3 == 0.0){
					summation = 0.0;
					for(t=0.0;t<=20.0;t+=1.0){
						if(l3-l2+t+m1<0.0 || l3-l1+t-m2<0.0 || l1+l2-l3-t<0.0 || l1-t-m1<0.0 || l2-t+m2 <0.0) continue;
						summation += pow(-1.0,t)/(SumPi(t)*SumPi(l3-l2+t+m1)*SumPi(l3-l1+t-m2)*SumPi(l1+l2-l3-t)*\
						SumPi(l1-t-m1)*SumPi(l2-t+m2));
					}
					diagonal = pow(-1.0,(double)(l1-l2-m3)) *\
					sqrt(SumPi(l1+l2-l3)*SumPi(l1-l2+l3)*SumPi(-l1+l2+l3)/SumPi(l1+l2+l3+1.0)) *\
					sqrt(SumPi(l1+m1)*SumPi(l1-m1)*SumPi(l2+m2)*SumPi(l2-m2)*SumPi(l3+m3)*SumPi(l3-m3)) *\
					summation;
					F6[(int)(m1+l1)][(int)(m2+l2)][(int)(m3+l3)] = diagonal;
}}}}}
