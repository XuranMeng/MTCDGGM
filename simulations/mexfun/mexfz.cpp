#include <Rcpp.h>
#include <math.h>
#include <vector>

// 定义 SQR 宏
#define SQR(x) ((x)*(x))

// normfun 函数
double normfun(const std::vector<double>& z, int kstart, int kend) {
    double nrm2 = 0.0;
    
    for (int i = kstart - 1; i < kend; ++i) {
        nrm2 += SQR(z[i]);
    }
    
    return sqrt(nrm2);
}

// [[Rcpp::export]]
double mexfz(const std::vector<double>& z, const std::vector<double>& ind, int grpNUM) {
    double fz = 0.0;
    
    for (int j = 0; j < grpNUM; ++j) {
        int kstart = (int) ind[3 * j];
        int kend = (int) ind[1 + 3 * j];
        double nrm = normfun(z, kstart, kend);
        fz += ind[2 + 3 * j] * nrm;
    }
    
    return fz;
}