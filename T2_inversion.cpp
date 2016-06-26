//T2 inversion functions
#include <complex>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <armadillo>
#include "NNLS.h"
#include "base_funcs.h"

//function reads in the real and imaginary data, and rotates them to output signal into R and noise into I.
void rotate_data(double * R, double * I, const int size)
{

double phi;
double sumR=0, sumI=0;
std::complex<double> m(0,-1);

for(unsigned int h=0; h<30;h++)
{
sumR=sumR+R[h];
sumI=sumI+I[h];
}
phi = atan2(sumI,sumR);

std::complex<double> Phi(phi,0);

for(unsigned int ii=0;ii<size;ii++)
{
std::complex<double> M(R[ii],I[ii]);
M=M*exp(m*Phi);
R[ii]=real(M);
I[ii]=imag(M);
}


}


//function "zeroes" the data by shifting it so its average is zero for the last 200 points. 
//Assumes the signal decays completely in the record.
void zero_data(double * R, double * I, const int size)
{
double m1=0,m2=0;
double a,b;
int count=0;

for(int N=size-200;N<size;N++)
{
m1=m1+R[N];
m2=m2+I[N];
count++;
}
a=m1/count;
b=m2/count;

for(unsigned int n=0;n<size;n++)
{
R[n]=R[n]-a;
I[n]=I[n]-b;
}
}

//Calculates the SNR of the dataset as (first data point)/(stdev of noise)
double SNR(double * R, double * I,const int size, double& stdev)
{
double var=0;
double m1=0;
double a;
int count=0,cnt=0;
double snr;

for(int N=size/2;N<size;N++)
{
m1=m1+R[N];
count++;
}
a=m1/count;

for(unsigned int n=size/2;n<size;n++)
{
var=var+(a-I[n])*(a-I[n]);
cnt++;
}
stdev=sqrt(var/cnt);
snr=R[1]/stdev;
return snr;
}

//Fills the vector T2 with n logarithmically spaced values from a to b
void logspace(double * T2, const double a, const double b, const int n)
{
double step=(b-a)/n;
double L;

for(unsigned int i=0;i<n;i++)
{
L=a+step*i;
T2[i]=pow(10,L);
}
}

//function to discard data after signal has decayed
void trunc_data(const double * R, const double * I, unsigned int& size)
{
	float m1,a,m2,b;
	int count=0,X=size/4;
	
	for(int N=size-1000;N<size;N++)
		{
		m1=m1+I[N];
		count++;
		}
	a=m1/count;
	count=0;
	
	do
	{
	for(int N=X;N<X+1000;N++)
		{
		m2=m2+R[N];
		count++;
		}
	b=m2/count;
	count=0;
	X=X+1000;
	m2=0;
	}
	while((std::abs(b)>std::abs(a*1.2))&&X<(size-1000));
	
	size=X+1000;
}

//function to cut the data 
void cut_data(const double * R, const double * I, const double * t, unsigned int& size, double * Rcut, double * Icut, double * tcut)
{
	int i, si=size/4, si2=size/8;
	
	for(i=0; i<si;i++) {Rcut[i]=R[i]; Icut[i]=I[i];tcut[i]=t[i];}
	for(i=0;i<si2;i++){Rcut[i+si]=R[si+2*i]; Icut[i+si]=I[si+2*i]; tcut[i+si]=t[si+2*i];}
	for(i=0;i<si2;i++){Rcut[si+si2+i]=R[si+si2+4*i]; Icut[si+si2+4*i]=I[si+si2+4*i]; tcut[si+si2+i]=t[si+si2+i];}
}
 

//T2 inversion based on whittall mod t2 fit code
void T2inversion(const double * T2, const double * R, const double * I, const double * time, const double alpha, double * T2dist, double * syn_data, const int size, const double stdev, const int nT2)
{
	const int k=nT2+1;
	int aa=((size+(k-3))*k);
	int bb=size*k;
	double* L=new double[aa];
	double* L0=new double[bb];
	
 	for(unsigned int i=0; i<size;i++)
 		{
 			for(unsigned int j=0; j<nT2;j++)
 				{
				//using row-wise mapping index (i,j) to nlines x k matrix
 				L[j+i*k]=exp(-time[i]/T2[j])/stdev;
 				L0[j+i*k]=exp(-time[i]/T2[j]);
 				}
 			L[k+i*k]=1/stdev;
 			L0[k+i*k]=1;
 		}
 
    	
 	for(unsigned int ii=size; ii<(size+(k-3));ii++)
 		{
 		for(unsigned int jj=0; jj<nT2;jj++)
 			{
 			if((ii-size)==jj) {L[jj+ii*k]=1*alpha;}
 			else if ((ii-size)==jj-1){L[jj+ii*k]=-2*alpha;}
 			else if ((ii-size)==jj-2){L[jj+ii*k]=1*alpha;}
 			else {L[jj+ii*k]=0;}
 			}
 		}
			
// m=NNLS(L,R); perform non-negative least-squares on L and R
//	std::cout<<(size+k-3)*k<<std::endl;
	NNLS(L,R,size,nT2,T2dist);

	//syn_data=L0*m
	matrix_mult(syn_data,L0,T2dist,size,k,k,1);
//end of function
}

