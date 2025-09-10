/***********************************************************************
* Compute Pz = ProxL2(z,c1,ind,grpNUM)
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
    double   *c1, *ind, *z, *Pz;
    double   cw, nrm; 

    ptrdiff_t  grpNUM;
    int   m, j, k,kstart, kend; 

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs > 4){
      mexErrMsgTxt("mexFnorm: requires at most 4 input arguments."); }
   if (nlhs > 1){ 
      mexErrMsgTxt("mexFnorm: requires at most 1 output argument."); }   

/* CHECK THE DIMENSIONS */
    
    
    m = mxGetM(prhs[0]); 
    z = mxGetPr(prhs[0]);
    c1 = mxGetPr(prhs[1]);
    ind = mxGetPr(prhs[2]);
    grpNUM = (int)*mxGetPr(prhs[3]);
    
    /***** create return argument *****/    
    plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL); 
    Pz = mxGetPr(plhs[0]);  

    /***** Do the computations in a subroutine *****/  
    for (j=0; j<grpNUM;j++) {
        cw = c1[0]*ind[2+3*j];
        kstart = (ptrdiff_t) ind[3*j];
        kend = (ptrdiff_t) ind[1+3*j];
        nrm = normfun(z, kstart, kend);
        /*if (j==0) { printf("1 nrm = %3.1f", nrm); }*/
        for (k=kstart-1; k< kend; k++) {
            if (nrm > cw) { 
               Pz[k] = z[k]*(1 - cw/nrm);
            } else{
               Pz[k] = 0.0;
            }
        }
    }
    return;
}
/**********************************************************/
// 输入参数：

// 	1.	z：输入向量，长度为 m，是需要投影的向量。
// 	2.	c1：标量，用于缩放组的L2范数约束。
// 	3.	ind：长度为 3 * grpNUM 的数组，描述每一组的索引和权重：
// 	•	ind[3 * j]：第 j 组的起始索引 kstart。
// 	•	ind[1 + 3 * j]：第 j 组的结束索引 kend。
// 	•	ind[2 + 3 * j]：第 j 组的权重 w_j。
// 	4.	grpNUM：整数，表示组的数量。

// 输出参数：

// 	1.	Pz：投影后的向量，长度为 m。

// 代码工作流程：

// 	1.	获取输入参数：
// 	•	z 是需要投影的输入向量。
// 	•	c1 是用于调整投影阈值的缩放因子。
// 	•	ind 包含每一组的起始和结束索引，以及权重。
// 	•	grpNUM 是组的数量。
// 	2.	初始化输出：
// 	•	创建输出向量 Pz 用于存储投影结果，初始时为零向量。
// 	3.	遍历每一组并计算L2范数：
// 	•	对每组 z[kstart:kend]，通过 normfun 函数计算L2范数。
// 	•	根据L2范数和阈值 cw = c1 * ind[2 + 3 * j] 比较：
// 	•	如果L2范数大于阈值，则将该组的向量 z[kstart:kend] 缩放到L2球内。
// 	•	如果L2范数小于或等于阈值，则将该组的所有元素置为零。
// 	4.	返回结果：
// 	•	最终输出的是投影后的向量 Pz。

// 假设输入如下：

// 	•	向量 z = [3, 4, 1, 2, 5, 6]，长度为 6。
// 	•	ind 描述了 2 个组：
// 	•	第 1 组：z[1:3] = [3, 4, 1]，权重 0.5。
// 	•	第 2 组：z[4:6] = [2, 5, 6]，权重 2.0。
// 	•	c1 = 0.75，用于调整阈值。
// 	•	grpNUM = 2，表示有两个组。

// 计算过程：

// 第 1 组（z[1:3] = [3, 4, 1]，权重 0.5）：

// 	•	计算L2范数：sqrt(3^2 + 4^2 + 1^2) = sqrt(9 + 16 + 1) = sqrt(26) ≈ 5.099。
// 	•	投影阈值：cw = 0.75 * 0.5 = 0.375。
// 	•	L2范数 5.099 大于 0.375，因此需要投影到L2球上：
// 	•	投影后的结果是：[3 * (1 - 0.375 / 5.099), 4 * (1 - 0.375 / 5.099), 1 * (1 - 0.375 / 5.099)] ≈ [2.779, 3.705, 0.926]。

// 第 2 组（z[4:6] = [2, 5, 6]，权重 2.0）：

// 	•	计算L2范数：sqrt(2^2 + 5^2 + 6^2) = sqrt(4 + 25 + 36) = sqrt(65) ≈ 8.062。
// 	•	投影阈值：cw = 0.75 * 2.0 = 1.5。
// 	•	L2范数 8.062 大于 1.5，需要投影：
// 	•	投影后的结果是：[2 * (1 - 1.5 / 8.062), 5 * (1 - 1.5 / 8.062), 6 * (1 - 1.5 / 8.062)] ≈ [1.627, 4.068, 4.881]。
