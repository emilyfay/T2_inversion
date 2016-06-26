# T2_inversion
C++ code to invert NMR decay data to calculate the distribution of T2 relaxation times

includes a makefile

run CPMG to calculate the T2 distribution from NMR relaxation data

Reads in data from a file, applies some simple pre-processing, then calculated a distribution of best-fit relaxation times using a non-negative least-squares algorithm with Tikohnov regularization. 
Inversion based on methodology described in:
Whittall, K. P., & MacKay, A. L. (1989). Quantitative interpretation of NMR relaxation data. Journal of Magnetic Resonance (1969), 84(1), 134-152.

*CPMG.cpp is currently written to load in the example data file C45.txt, but includes code to read in a user-specified data file
The data file should contain entries for time (in seconds) in the first column, the real signal in the second column, and the imaginary signal in the third column (in arbitrary units)
