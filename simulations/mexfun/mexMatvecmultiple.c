#include "mex.h"

/* 
 * MEX function to compute fun = m' * Sigma_hat * m
 * Input:
 *   - Sigma_hat: the sample covariance matrix (p x p)
 *   - m: the vector (p x 1)
 * Output:
 *   - result: the scalar result of m' * Sigma_hat * m
 */
void compute_function(double* Sigma_hat, double* m, int p, double* result) {
    double sum = 0.0;

    // Perform the matrix-vector multiplication m' * Sigma_hat * m
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            sum += m[i] * Sigma_hat[i * p + j] * m[j];
        }
    }

    *result = sum;
}

/* The gateway function for MATLAB to call */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "Two inputs required.");
    }
    
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "One output required.");
    }
    
    /* Get the inputs */
    double *Sigma_hat = mxGetPr(prhs[0]);
    double *m = mxGetPr(prhs[1]);
    
    /* Get the dimensions */
    int p = mxGetN(prhs[0]);
    
    /* Create the output variable */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *result = mxGetPr(plhs[0]);
    
    /* Call the computational function */
    compute_function(Sigma_hat, m, p, result);
}