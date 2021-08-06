//split data used//
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

int main(int argc,char* argv[]){
	int i,Data;
	double col[12],max=-20.0,min=20.0;
	FILE *FpR,*FpW;
	if(argc!=4){
		printf("argv should be 4");
		exit (1);
	}
	FpR = fopen(argv[1],"r");
	FpW = fopen(argv[2],"a");	
	Data = atoi(argv[3]);
	while(fscanf(FpR,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",\
	&col[0],&col[1],&col[2],&col[3],&col[4],&col[5],&col[6],&col[7],&col[8],&col[9],&col[10],&col[11]) != EOF){
		if(col[Data]>max) max = col[Data];
		if(col[Data]<min) min = col[Data];
	}
	fprintf(FpW,"%s %d max is %lf min is %lf\n",argv[1],Data,max,min);
	fclose(FpR);
	fclose(FpW);
	return 0;
}
