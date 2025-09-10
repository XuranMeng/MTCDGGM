/***********************************************************************
* Compute [Pz,indicate] = ProjL2(z,c1,ind,grpNUM)
* by Xudong Li, 8 May 2017
***********************************************************************/

#include <math.h>
#include <mex.h>
#include <matrix.h>
#include <string.h> /* needed for memcpy() */


#if !defined(SQR)
#define SQR(x) ((x)*(x))
#endif


/**********************************************************
* 
***********************************************************/
double normfun(const double *z, int kstart, int kend)

{ 
    double nrm, nrm2 =0.0;
    ptrdiff_t i;
    
    for (i=kstart-1; i< kend; i++) {
        nrm2 += SQR(z[i]);
    }
    nrm = sqrt(nrm2);
    /*printf(" kstart = %d, kend = %d ", kstart, kend);*/
    return nrm;
}


void mexFunction(
      int nlhs,   mxArray  *plhs[], 
      int nrhs,   const mxArray  *prhs[] )

{        
    double   *c1, *ind, *z, *Pz, *indicate;
    double   cw, nrm;

    ptrdiff_t  grpNUM;
    int   m, j, k,kstart, kend; 

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs > 4){
      mexErrMsgTxt("mexFnorm: requires at most 4 input arguments."); }
   if (nlhs > 2){ 
      mexErrMsgTxt("mexFnorm: requires at most 2 output argument."); }   

/* CHECK THE DIMENSIONS */
    
    
    m = mxGetM(prhs[0]); 
    z = mxGetPr(prhs[0]);
    c1 = mxGetPr(prhs[1]);
    ind = mxGetPr(prhs[2]);
    grpNUM = (int)*mxGetPr(prhs[3]);
    
    /***** create return argument *****/    
    plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(grpNUM,1,mxREAL);
    Pz = mxGetPr(plhs[0]);  
    indicate = mxGetPr(plhs[1]); 

    /***** Do the computations in a subroutine *****/  
    for (j=0; j<grpNUM;j++) {
        cw = c1[0]*ind[2+3*j];
        kstart = (ptrdiff_t) ind[3*j];
        kend = (ptrdiff_t) ind[1+3*j];
        nrm = normfun(z, kstart, kend);
        /*if (j==0) { printf("1 nrm = %3.1f", nrm); }*/
        if (nrm > cw) {
            indicate[j] = nrm;
        } else{
            indicate[j] = 0;
        }
        for (k=kstart-1; k< kend; k++) {
            if (nrm > cw) { 
               Pz[k] = z[k]*(cw/nrm);
            } else{
               Pz[k] = z[k];
            }
        }
    }
    return;
}
/**********************************************************/


// 输入参数：

// 	1.	z：向量，长度为 m，这是需要进行投影的向量。
// 	2.	c1：标量，定义投影时的L2球的缩放因子。
// 	3.	ind：描述组结构的矩阵，维度为 (3 * grpNUM)：
// 	•	ind[3*j]：每一组的起始索引 kstart。
// 	•	ind[1 + 3*j]：每一组的结束索引 kend。
// 	•	ind[2 + 3*j]：每一组的权重 w_j。
// 	4.	grpNUM：整数，表示有多少组。

// 输出参数：

// 	1.	Pz：投影后的向量，长度为 m。
// 	2.	indicate：一个向量，长度为 grpNUM，存储每一组的L2范数。如果该组的L2范数超过了设定的阈值 cw，那么返回范数，否则返回0。


// 	•	向量 z = [2, 3, 0.1, 0.2, 0.3, 5, 6, 1]。
// 	•	其他参数不变。

// 计算过程：

// 第 1 组：

// 同之前一样，结果为 [0.208, 0.312]，indicate[1] = 3.605。

// 第 2 组（z[3:5] = [0.1, 0.2, 0.3]，权重 2.0）：

// 	•	计算L2范数：sqrt(0.1^2 + 0.2^2 + 0.3^2) = sqrt(0.01 + 0.04 + 0.09) = sqrt(0.14) ≈ 0.374。
// 	•	投影阈值：cw = 1.5。
// 	•	L2范数 0.374 小于 1.5，不需要投影。
// 	•	投影结果为 Pz[3:5] = [0.1, 0.2, 0.3]，indicate[2] = 0。

// 第 3 组：

// 同之前一样，结果为 [0.477, 0.573, 0.095]，indicate[3] = 7.874。

