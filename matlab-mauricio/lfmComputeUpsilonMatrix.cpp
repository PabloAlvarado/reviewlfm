/*
 * lfmComputeUpsilonMatrix.cpp 
 * 
 * Computes Upsilon Matrix.
 * inputs:
 *  gamma: gamma value system
 *  sigma2: squared lengthscale
 *  t1: first time input (x 1)
 *  t2: second time input (y 1)
 * outputs:
 *  Upsilon Matrix (x y)
 *
 * This is a MEX-file for MATLAB.
 *
 * Diego Agudelo, 2015.
 *
 * How to compile this?
 *
 * 1) Generate the object file for Faddeeva.cc:
 *  > mex -c Faddeva.cc
 * 2) Compile the mex c++ file and link with the Faddeeva library
 *  > mex lfmComputeUpsilonMatrix.cpp Faddeeva.o
*/

#include <iostream>
#include "Faddeeva.hh"
#include "mex.h"

using namespace std;

bool isScalar(const mxArray *element) {
  mwSize rows = mxGetM(element);
  mwSize cols = mxGetN(element);
  return (rows == 1 && cols == 1);
}

bool isColumnVector(const mxArray *element) {
  mwSize cols = mxGetN(element);
  return (cols == 1);
}

#define MAX_COLS 10000

void computeUpsilon(
    complex<double> gamma,
    double sigma2,
    double *t1,
    double *t2,
    mwSize rows,
    mwSize cols,
    double *upsilonReal,
    double *upsilonComplex) {
  complex<double> complexJ(0, 1);
  complex<double> complexMinusJ(0, -1);
  double dev = sqrt(sigma2);

  // assert(cols <= MAX_COLS);
  complex<double> Z2[MAX_COLS];
  complex<double> WOFZ2[MAX_COLS];
  
  // complex<double> Z2[cols];
  // complex<double> WOFZ2[cols];
  
  for(mwSize j = 0; j < cols; j++) {
    Z2[j] = t2[j] / dev + dev * gamma / 2.0;
    if (real(Z2[j]) >= 0.0)
      WOFZ2[j] = Faddeeva::w(complexJ * Z2[j]);
    else
      WOFZ2[j] = Faddeeva::w(complexMinusJ * Z2[j]);
  }
  
  for(mwSize i = 0; i < rows; i++) {
    for(mwSize j = 0; j < cols; j++) {
      double t1MinusT2 = t1[i] - t2[j];
      complex<double> Z1 = t1MinusT2 / dev - dev * gamma / 2.0;
      complex<double> WOFZ1;
      if (real(Z1) >= 0.0)
        WOFZ1 = Faddeeva::w(complexJ * Z1);
      else
        WOFZ1 = Faddeeva::w(complexMinusJ * Z1);
      complex<double> ans;
      if (real(Z1) >= 0.0 && real(Z2[j]) >= 0.0) {
        ans = (
          2.0 * exp(sigma2 * (gamma*gamma)/4.0 - gamma * t1MinusT2)
          - exp( -(t1MinusT2 * t1MinusT2 / sigma2) + log(WOFZ1))
          - exp( -(t2[j] * t2[j] / sigma2) - gamma * t1[i] + 
          log(WOFZ2[j]) )
        ); 
      }
      if (real(Z1) < 0.0 && real(Z2[j]) >= 0.0) {
        ans = (
          exp( -(t1MinusT2 * t1MinusT2 / sigma2) + log(WOFZ1))
          - exp( -(t2[j] * t2[j] / sigma2) - gamma * t1[i] + 
          log(WOFZ2[j]) )
        );  
      }
      if (real(Z1) >= 0.0 && real(Z2[j]) < 0.0) {
        ans = (
          - exp( -(t1MinusT2 * t1MinusT2 / sigma2) + log(WOFZ1))
          + exp( -(t2[j] * t2[j] / sigma2) - gamma * t1[i] + 
          log(WOFZ2[j]) )
        );
      }
      if (real(Z1) < 0.0 && real(Z2[j]) < 0.0) {
        ans = (
          -2.0 * exp(sigma2 * (gamma*gamma)/4.0 - gamma * t1MinusT2)
          + exp( -(t1MinusT2 * t1MinusT2 / sigma2) + log(WOFZ1))
          + exp( -(t2[j] * t2[j] / sigma2) - gamma * t1[i] + 
          log(WOFZ2[j]) )
        );
      }
      upsilonReal[i + rows * j] = real(ans);
      upsilonComplex[i + rows * j] = imag(ans);
    }
  }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Input validation
  if (nrhs != 4)
    mexErrMsgTxt("Four inputs required.");
  if (!isScalar(prhs[0]))
    mexErrMsgTxt("Input #1 is not a scalar.");
  if (!isScalar(prhs[1]))
    mexErrMsgTxt("Input #2 is not a scalar.");
  if (!isColumnVector(prhs[2]))
    mexErrMsgTxt("Input #3 is not a column vector.");
  if (!isColumnVector(prhs[3]))
    mexErrMsgTxt("Input #4 is not a column vector.");
  // Converting MATLAB variables to C++ variables
  double *gammaReal = mxGetPr(prhs[0]);
  double *gammaComplex = mxGetPi(prhs[0]);
  complex<double> gamma = gammaReal[0];
  if (mxIsComplex(prhs[0]))
    gamma = complex<double>(gammaReal[0], gammaComplex[0]);
  double sigma2 = mxGetScalar(prhs[1]);
  double *t1 = mxGetPr(prhs[2]);
  double *t2 = mxGetPr(prhs[3]);
  mwSize rows = mxGetM(prhs[2]);
  mwSize cols = mxGetM(prhs[3]);
  // Creating output matrix
  plhs[0] = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
  double *upsilonReal = mxGetPr(plhs[0]);
  double *upsilonComplex = mxGetPi(plhs[0]);
  computeUpsilon(gamma, sigma2, t1, t2, rows, cols, upsilonReal, upsilonComplex);
}
