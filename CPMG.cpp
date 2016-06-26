//read in CPMG data from the Maran

#include <iostream>
#include <stdio.h>
#include <string>
#include "T2_inversion.h"
#include <math.h>
#include <armadillo>

using namespace std;

int main()
{
/* uncomment section below to get filename from user */
	//cout<<"Input CPMG data filename:"<<endl;
	//string filename;
	//cin >> filename;
	//const char * fn = filename.c_str();
const char * fn ="C45.txt"; //comment out if using user-input filename
/*load file*/ 
FILE* f=fopen(fn,"r"); //choose an output file from a CPMG experiment on the Maran
if(f==NULL)
{printf("There was an error loading this file.\n");
exit(EXIT_FAILURE);}

unsigned int nlines=0;
int ch;

	while(EOF!=(ch=getc(f)))
	if('\n'==ch)
	++nlines;
	
	rewind(f);
	printf("%i\n",nlines);

int ni=0, nt=0, nR=0, nI=0; //initialize counters
double* time=new double[nlines]; 
double* R=new double[nlines];
double* I=new double[nlines];

	while(!feof(f)&&ni<(3*nlines)) //until the end of file is reached
	{
	
	if(ni%3==0) //sort first column into time vector
	{
	fscanf(f,"%lf", &time[nt]); 
	nt++;
	}
	else if(ni%3==1) //sort second column into Real data vector
	{
	fscanf(f,"%lf", &R[nR]); 
	nR++;
	}
	else if(ni%3==2) //sort third column into Imaginary data vector
	{
	fscanf(f,"%lf", &I[nI]); 
	nI++;
	}
	
	ni++;
	}
	printf("Read in %i time points, %i real data points and %i imaginary data points.\n",nt,nR,nI);
	fclose(f);

/* apply functions to rotate the data to maximize the real component and adjust for any offset */
rotate_data(R,I,nlines);
zero_data(R,I,nlines);

/* calculate the SNR of the data */
double stdev=0;
double snr=SNR(R,I,nlines,stdev);
printf("SNR = %f \n stdev=%f \n",snr,stdev);

/* truncate the data to speed up the inversion */
trunc_data(R,I,nlines);

unsigned int scut=nlines/4+nlines/8+nlines/8;
double tcut[scut], Rcut[scut],Icut[scut];

cut_data(R,I,time,nlines,Rcut,Icut,tcut);

delete[] R; delete[] I; delete[] time;

int nT2=200;
double T2[nT2];
logspace(T2,-4,0.7,nT2);

double alpha = 1;
double  T2dist[nT2+1];
double syn_data[scut];

/* calculate the relaxation time distribution using a NNLS algorithm */
printf("Beginning inversion");
T2inversion(T2,Rcut,Icut,tcut,alpha,T2dist,syn_data,scut,stdev,nT2);

/* Save output if requested by user*/
cout<<"Save output? Type Y to save, anything else to exit without saving."<<endl;
char Save[1];
cin >> Save;

if(strcmp(Save,"Y")==0||strcmp(Save,"y")==0)
{
char outfile[50];
sprintf(outfile,"Data_%s",fn);
FILE* fout=fopen(outfile,"wb");
for(unsigned i=0;i<scut;i++)
fprintf(fout,"%f\t%f\t%f\n",tcut[i],Rcut[i],Icut[i]);
fclose(fout);

char outfileT2[50];
sprintf(outfileT2,"T2dist_%s",fn);
FILE* fout2=fopen(outfileT2,"wb");
for(unsigned i=0;i<(nT2+1);i++)
fprintf(fout2,"%f\n",T2dist[i]);
fclose(fout2);

char outfileSyn[50];
sprintf(outfileSyn,"Syn_data_%s",fn);
FILE* fout3=fopen(outfileSyn,"wb");
for(unsigned i=0;i<(scut);i++)
fprintf(fout3,"%f\t%f\n",tcut[i],syn_data[i]);
fclose(fout3);


printf("Saved \n");
}
return 0;
}