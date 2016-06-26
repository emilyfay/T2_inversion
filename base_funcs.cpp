#include <stdio.h>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//function for multiplying two matrices
void matrix_mult(double * product, const double * matA, const double * matB, const int m1, const int n1, const int m2, const int n2)
{
if(n1!=m2)
{ printf("Error, cannot multiply these matrices");
exit(EXIT_FAILURE);
}
double sp=0;
	for(unsigned int i=0;i<m1;i++)
	{
		for(unsigned int j=0; j<n2;j++)
		{
			for(unsigned int jj=0; jj<n1;jj++)
			{
			sp=matA[i*n1+jj]*matB[jj*n2+j]+sp;
			}
		product[i*n2+j]=sp;
		sp=0;
		}
	}


}
////////////////////////////////////////////////////////////////////////////////
//function for transposing a matrix
void matrix_trans(double * AT, const double * A, const int m, const int n)
{
	for(unsigned int i=0; i<m;i++)
	{
		for(unsigned int j=0; j<n;j++)
		{
			AT[j*m+i]=A[i*n+j];
		}
	}

}
////////////////////////////////////////////////////////////////////////////////
//function to determine if any values in a vector are true
bool any_true(const bool * IN, const int size)
{
bool B=false;
	for(unsigned int i=0;i<size;i++)
	{
		if(IN[i])
		{B=true; break;}
	}
return B;
}

////////////////////////////////////////////////////////////////////////////////
//function to determine if values in a vector determined by the indices are greater than the test
bool any_greater(const double * IN, const bool * indices, const int size, const double test)
{
bool B=false;

	for(unsigned int i=0;i<size;i++)
	{
		if(indices[i]&&(IN[i]>test))
		{B=true; break;}
	}
return B;
}
////////////////////////////////////////////////////////////////////////////////
//function to fill a vector of length n with chosen value
void fill(double * X, const int size, const double x)
{
	for(unsigned int i=0;i<size;i++)
	{
	X[i]=x;
	}
}
////////////////////////////////////////////////////////////////////////////////
//function to fill a vector of length n with chosen value
void fill(bool* X, const int size, const bool x)
{
	for(unsigned int i=0;i<size;i++)
	{
	X[i]=x;
	}
}
////////////////////////////////////////////////////////////////////////////////
//function to fill a vector of length n with chosen value
void fill(int* X, const int size, const int x)
{
	for(unsigned int i=0;i<size;i++)
	{
	X[i]=x;
	}
}

////////////////////////////////////////////////////////////////////////////////
//function to fill selected indices of a vector of length n with chosen value
void fillI(int* X, const bool* indices, const int size, const int x)
{
	for(unsigned int i=0;i<size;i++)
	{
	if(indices[i]) {X[i]=x;}
	}
}

////////////////////////////////////////////////////////////////////////////////
//function to fill selected indices of a vector of length n with chosen value
void fillI(double* X, const bool* indices, const int size, const double x)
{
	for(unsigned int i=0;i<size;i++)
	{
	if(indices[i]) {X[i]=x;}
	}
}

////////////////////////////////////////////////////////////////////////////////
//function to fill selected indices of a vector of length n with chosen value
void fillI(double* X, const bool* indices, const int size, const double* x)
{
	for(unsigned int i=0;i<size;i++)
	{
	if(indices[i]) {X[i]=x[i];}
	}
}

////////////////////////////////////////////////////////////////////////////////
// function to return the index of the maximum entry 
int indMAX(const double* X, const int size)
{
int imax=0;
double X0=X[0];
	for(unsigned int i=1;i<size;i++)
	{
	if(X[i]>X0){imax=i; X0=X[i];}
	}
	return imax;
}

