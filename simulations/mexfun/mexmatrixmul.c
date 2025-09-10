#include "mex.h"
#include <math.h>
#include <matrix.h>
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of input and output arguments */
    if(nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidNumInputs", 
                          "Two input matrices are required.");
    }
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidNumOutputs", 
                          "One output matrix is required.");
    }

    /* Get the dimensions of the input matrices A and B */
    mwSize m = mxGetM(prhs[0]);   // Number of rows of A
    mwSize k = mxGetN(prhs[0]);   // Number of columns of A (also rows of B)
    mwSize n = mxGetN(prhs[1]);   // Number of columns of B

    /* Check if the inner dimensions match for matrix multiplication */
    if(mxGetM(prhs[1]) != k) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:innerDimMismatch", 
                          "Inner dimensions of matrices do not match.");
    }

    /* Create the output matrix C (m x n) */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);

    /* Get pointers to the data in A, B, and C */
    double *A = mxGetPr(prhs[0]);
    double *B = mxGetPr(prhs[1]);
    double *C = mxGetPr(plhs[0]);

    /* Perform matrix multiplication C = A * B */
    for(mwSize i = 0; i < m; i++) {
        for(mwSize j = 0; j < n; j++) {
            C[i + j*m] = 0;  // Initialize C[i,j]
            for(mwSize l = 0; l < k; l++) {
                C[i + j*m] += A[i + l*m] * B[l + j*k];
            }
        }
    }
}