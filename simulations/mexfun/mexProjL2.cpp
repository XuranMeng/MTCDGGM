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
Rcpp::List mexProjL2(const std::vector<double>& z, double c1, const std::vector<double>& ind, int grpNUM) {
    int m = z.size();
    std::vector<double> Pz(m, 0.0);
    std::vector<double> indicate(grpNUM, 0.0);
    
    for (int j = 0; j < grpNUM; ++j) {
        double cw = c1 * ind[2 + 3 * j];
        int kstart = (int) ind[3 * j];
        int kend = (int) ind[1 + 3 * j];
        double nrm = normfun(z, kstart, kend);

        // 判断是否需要进行投影
        if (nrm > cw) {
            indicate[j] = nrm;
        } else {
            indicate[j] = 0.0;
        }

        // 计算 Pz 向量
        for (int k = kstart - 1; k < kend; ++k) {
            if (nrm > cw) {
                Pz[k] = z[k] * (cw / nrm);
            } else {
                Pz[k] = z[k];
            }
        }
    }

    return Rcpp::List::create(Rcpp::Named("Pz") = Pz, 
                              Rcpp::Named("indicate") = indicate);
}