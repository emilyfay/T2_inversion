//header template file

#ifndef T2_INVERSION_H
#define T2_INVERSION_H


//declarations
void rotate_data(double * R, double * I, const int size);
void zero_data(double * R, double * I, const int size);
double SNR(double * R, double * I,const int size, double& stdev);
void logspace(double * T2, const double a, const double b, const int n);
void trunc_data(const double * R, const double * I, unsigned int& size);
void cut_data(const double * R, const double * I, const double * t, unsigned int& size, double * Rcut, double * Icut, double * tcut);
void T2inversion(const double * T2, const double * R, const double * I, const double * time, const double alpha, double * T2dist, double * syn_data, const int size, const double stdev, const int nT2);

#endif