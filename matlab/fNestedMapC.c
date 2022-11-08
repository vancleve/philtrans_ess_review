/*
 * =============================================================
 * MEX file
 *
 * Jeremy Van Cleve <vancleve@stanford.edu>
 * =============================================================
 */
 
#include "mex.h"
#include <math.h>

double fMapC(double x, double mu, double w1, double w3)
{
  return ((1-mu)*w1*x + mu*w3*(1-x)) / (w1*x + w3*(1-x));
}

void fNestedMapC(double xin, double *xout, double mu, double w11, double w31, double w12, double w32, int n)
{
  int i;
  
  *xout = xin;
  for (i=0; i<n; i++) {
    *xout = fMapC(*xout, mu, w11, w31);
  }
  for (i=0; i<n; i++) {
    *xout = fMapC(*xout, mu, w32, w12);
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double xin, mu, w11, w31, w12, w32, n;
  double *xout;

  /* Check for correct number of input and output arguments */
  if (nrhs < 7) {
    mexErrMsgTxt("This function requires five arguments.");
  }
  if (nlhs != 1) {
    mexErrMsgTxt("This function requires one output.");
  }
  
  /* Assign each input and output. */
  xin = mxGetScalar(prhs[0]);
  mu = mxGetScalar(prhs[1]);
  w11 = mxGetScalar(prhs[2]);
  w31 = mxGetScalar(prhs[3]);
  w12 = mxGetScalar(prhs[4]);
  w32 = mxGetScalar(prhs[5]);
  n = mxGetScalar(prhs[6]);

  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);

  /* Pointer to output */
  xout = mxGetPr(plhs[0]);

  /* Call the subroutine. */
  fNestedMapC(xin, xout, mu, w11, w31, w12, w32, n);
}
