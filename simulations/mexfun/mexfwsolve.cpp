#include <Rcpp.h>
#include <math.h>
#include <vector>

// [[Rcpp::export]]
std::vector<double> mexfwsolve(const Rcpp::NumericMatrix& R, const Rcpp::NumericVector& b) {
  int n = b.size();

  // 检查矩阵是否为方阵
  if (R.nrow() != n || R.ncol() != n) {
    Rcpp::stop("mexfwsolve: R should be square.");
  }

  // 检查是否为上三角矩阵（保留原有的判断逻辑）
  if (R(0, 0) == 0) {
    Rcpp::stop("mexfwsolve: R not upper triangular.");
  }

  // 创建结果向量
  std::vector<double> x(n, 0.0);

  // 求解第一个方程
  x[0] = b[0] / R(0, 0);

  // 前代法求解其余方程
  for (int j = 1; j < n; ++j) {
    double tmp = 0.0;
    for (int k = 0; k < j; ++k) {
      tmp += R(k, j) * x[k];
    }
    x[j] = (b[j] - tmp) / R(j, j);
  }

  return x;
}