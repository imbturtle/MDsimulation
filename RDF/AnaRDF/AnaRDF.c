#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(int argc,char* argv[]){
	int i;
	double col[2],max,min,rmin,rmax;
	FILE *fr,*fw;
	if(argc!=4){
		printf("argv should be 4");
		exit (1);
	}
	fr = fopen(argv[1],"r");
	fw = fopen(argv[2],"a");	
	fprintf(fw,"%s ",argv[3]);
	max=0.0;
	min=20.0;
	while(fscanf(fr,"%lf %lf",&col[0],&col[1]) != EOF){//find maximun
		if(col[1] > max) {
			rmax = col[0];
			max = col[1];
		}
		if(col[1] < max) break;
	}	
	fprintf(fw,"%.3lf %.3lf ",rmax,max);
	while(fscanf(fr,"%lf %lf",&col[0],&col[1]) != EOF){//find mimun
		if(col[1] < min) {
			rmin = col[0];
			min = col[1];
		}	
		if(col[1] > min) break;
	}
	fprintf(fw,"%.3lf %.3lf %.3lf ",rmin,min,min/max);
	max=0.0;
	min=20.0;
	while(fscanf(fr,"%lf %lf",&col[0],&col[1]) != EOF){//find maximun
		if(col[1] > max) {
			rmax = col[0];
			max = col[1];
		}
		if(col[1] < max) break;
	}
	fprintf(fw,"%.3lf %.3lf ",rmax,max);
	while(fscanf(fr,"%lf %lf",&col[0],&col[1]) != EOF){//find mimun
		if(col[1] < min) {
			rmin = col[0];
			min = col[1];
		}	
		if(col[1] > min) break;
	}
	fprintf(fw,"%.3lf %.3lf\n",rmin,min);
	fclose(fr);
	fclose(fw);
	return 0;
}
