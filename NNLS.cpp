#include <complex>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include "base_funcs.h"
#include <armadillo>

using namespace arma;

int indMax(vec X, const int size);
////////////////////////////////////////////////////////////////////////////////
// non-negative least-squares algorithm
void NNLS(const double * LL, const double * R, const int size, const int nT2, double * T2dist)
{
//size of R should be size x 1
//size of L should be (size+(nT2+1)-3) x (nT2+1)
const int n=nT2+1;
const int N=size+n-3;
mat x(n,1); x.zeros();
ivec P(n); P.zeros(); //positive set
ivec Z(n); Z.ones(); //zero set
vec d(N), z(n), wz(n), Q(n);
vec INF(n); INF.fill(-1*datum::inf);
mat L(N,n);
const double tol= 1e-4;//1.5e-8;
//set up iteration criterion and counters
int INit=0, OUTit=0, MAXit=3*n, tt, ii, i, j;
 
 //initialize data vector from R data with zeros added on at end
for( i=0;i<(size);i++)
{
	d[i]=R[i];
}
for(i=size;i<(N);i++)
{
	d[i]=0;
}

//convert array LL into matrix L
for(i=0; i<N;i++)
{
for(j=0;j<n;j++)
L(i,j)=LL[i*n+j];
}

//calulate initial value for residual vector
vec resid=d-(L*x);
vec w=L.t()*resid;

//iterate while there are still elements in the zero set and while those elements of the w matrix are larger than the tolerance
while(any(Z) && any((w%Z)>tol))
	{
	z.zeros(); // reset the intermediate solution
	//create wz, a Lagrange multiplier vector of variables in the zero set
	wz=P%INF; //equivalent to wz(P)=-Inf
	for(i=0;i<n;i++) //loop to set wz(Z)=w(Z);
	{
	if(Z(i)>0.5) wz(i)=w(i);
	}
	
	tt=indMax(wz,n); //returns the index of the variable with the largest Lagrange multiplier
	// move this variable from the zero set to the positive set
	P(tt)=1;
	Z(tt)=0;
	
	//calculate the intermediate solution using only variables in the positive set
	int Sx=sum(P);
	mat Lp(N,Sx);
	vec zp(Sx);
	ii=0;
	for(i=0;i<n;i++)
		{ if(P(i)>0.5) {Lp.col(ii)=L.col(i);ii++;} }

	zp =solve(Lp,d);
	ii=0;
	for(i=0;i<n;i++)
		{ if(P(i)>0.5) {z(i)=zp(ii);ii++;} }

//begin the inner loop to remove elements from the positive set which no longer belong
		while(any(zp<=0))
		{
			INit++;
			//stop and return z if the maximum numbr of iterations has been reached
			if(INit>MAXit) { for(ii=0;ii<n;ii++){T2dist[ii]=z(ii);} std::cout<<"Exceeded max iterations"<<std::endl; return;}
			
			//find the indices where z is approximately negative	
			for(uint jj=0; jj<n;jj++)
				{if((z(jj)<=0)&&(P(jj)>0.5))	Q(jj)=x(jj)/(x(jj)-z(jj));
				else Q(jj)=datum::inf;}
			//choose new x to keep new x nonnegative
			double alpha=min(Q);
			x = x+alpha*(z-x);
			//Reset Z and P given intermediate values of x
			for(uint jj=0; jj<n;jj++)
				{ if((std::abs(x(jj))<tol && P(jj)>0.5) || Z(jj)>0.5) {Z(jj)=1; P(jj)=0;}
				else {Z(jj)=0; P(jj)=1;}
				}
			z.zeros();
		
			//re-solve for z
			Sx=sum(P);
			Lp.reset(); Lp.set_size(N,Sx);
			zp.reset(); zp.set_size(Sx);
			ii=0;
			for(i=0;i<n;i++)
				{ if(P(i)>0.5) {Lp.col(ii)=L.col(i);ii++;} }
		
			zp =solve(Lp,d);
			ii=0;
			for(i=0;i<n;i++)
				{ if(P(i)>0.5) {z(i)=zp(ii);ii++;} }

		}
	x=z;
	resid = d-(L*x);
	w=L.t()*resid;

	OUTit++;
	}
//transfer z vector into T2dist array to return
for(ii=0;ii<n;ii++){T2dist[ii]=z(ii);}
std::cout<<"Met stopping criteria"<<std::endl;
//end of function
}

int indMax(vec X, const int size)
{
int imax=0;
double X0=X(0);
	for(unsigned int i=1;i<size;i++)
	{
	if(X(i)>X0){imax=i; X0=X(i);}
	}
	return imax;
}
