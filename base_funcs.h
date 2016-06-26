#ifndef BASE_FUNCS_H
#define BASE_FUNCS_H


void matrix_mult(double * product, const double * matA, const double * matB, const int m1, const int n1, const int m2, const int n2);
void matrix_trans(double * AT, const double * A, const int m, const int n);
bool any_true(const bool * IN, const int size);
bool any_greater(const double * IN, const bool * indices, const int size, const double test);
void fill(double * X, const int size, const double x);
void fill(bool* X, const int size, const bool x);
void fill(int* X, const int size, const int x);
void fillI(double* X, const bool* indices, const int size, const double x);
void fillI(int* X, const bool* indices, const int size, const int x);
void fillI(double* X, const bool* indices, const int size, const double* x);
int indMAX(const double* X, const int size);

#endif