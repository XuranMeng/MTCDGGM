#include <Rcpp.h>
#include <math.h>
#include <vector>

// 定义 SQR 宏
#define SQR(x) ((x)*(x))

// realdot 函数：计算两个向量的点积，使用循环展开优化
double realdot(const std::vector<double>& x, const std::vector<double>& y, const int n) {
    double r = 0.0;
    int i = 0;
    
    for (; i < n - 7; i += 8) {  // 8级展开
        r += x[i] * y[i];
        r += x[i+1] * y[i+1];
        r += x[i+2] * y[i+2];
        r += x[i+3] * y[i+3];
        r += x[i+4] * y[i+4];
        r += x[i+5] * y[i+5];
        r += x[i+6] * y[i+6];
        r += x[i+7] * y[i+7];
    }
    
    for (; i < n; ++i) {  // 处理剩余部分
        r += x[i] * y[i];
    }

    return r;
}

// subscalarmul 函数：x -= alpha * y，使用循环展开
void subscalarmul(std::vector<double>& x, const double alpha, const std::vector<double>& y, const int n) {
    int i = 0;
    
    for (; i < n - 7; i += 8) {  // 8级展开
        x[i] -= alpha * y[i];
        x[i+1] -= alpha * y[i+1];
        x[i+2] -= alpha * y[i+2];
        x[i+3] -= alpha * y[i+3];
        x[i+4] -= alpha * y[i+4];
        x[i+5] -= alpha * y[i+5];
        x[i+6] -= alpha * y[i+6];
        x[i+7] -= alpha * y[i+7];
    }
    
    for (; i < n; ++i) {  // 处理剩余部分
        x[i] -= alpha * y[i];
    }
}

// lbsolve 函数：解 y，从 U' * y = x 得出
void lbsolve(std::vector<double>& y, const std::vector<double>& U, const std::vector<double>& x, const int n) {
    for (int k = 0; k < n; ++k) {
        y[k] = (x[k] - realdot(y, std::vector<double>(U.begin() + k * n, U.begin() + (k + 1) * n), k)) / U[k * n + k];
    }
}

// ubsolve 函数：解 xnew，从 U * xnew = x 得出
void ubsolve(std::vector<double>& x, const std::vector<double>& U, const int n) {
    for (int j = n - 1; j >= 0; --j) {
        x[j] /= U[j * n + j];
        subscalarmul(x, x[j], std::vector<double>(U.begin() + j * n, U.begin() + (j + 1) * n), j);
    }
}

// [[Rcpp::export]]
std::vector<double> mextriang(const Rcpp::NumericMatrix& U, const Rcpp::NumericVector& b, int options = 1) {
    int n = U.nrow();
    
    if (U.nrow() != U.ncol()) {
        Rcpp::stop("U must be a square matrix.");
    }
    
    if (b.size() != n) {
        Rcpp::stop("Size of U and b mismatch.");
    }

    std::vector<double> y(n);
    std::vector<double> Uvec(U.begin(), U.end());
    std::vector<double> bvec(b.begin(), b.end());

    // 选择执行 ubsolve 或 lbsolve
    if (options == 1) {
        y = bvec;
        ubsolve(y, Uvec, n);
    } else if (options == 2) {
        lbsolve(y, Uvec, bvec, n);
    } else {
        Rcpp::stop("Invalid option for solving.");
    }
    
    return y;
}