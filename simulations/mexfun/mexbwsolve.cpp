#include <Rcpp.h>
#include <math.h>
#include <vector>

// [[Rcpp::export]]
std::vector<double> mexbwsolve(const Rcpp::NumericMatrix& Rt, const Rcpp::NumericVector& b) {
  int n = b.size();

  // 检查矩阵维度是否为方阵
  if (Rt.nrow() != n || Rt.ncol() != n) {
    Rcpp::stop("mexbwsolve: Rt should be square.");
  }

  // 检查是否为下三角矩阵（保留之前代码的逻辑）
  if (Rt(n-1, n-1) == 0) {
    Rcpp::stop("mexbwsolve: Rt is not lower triangular.");
  }

  // 创建结果向量
  std::vector<double> x(n, 0.0);

  // 最后一个方程：x[n-1] = b[n-1] / Rt[n-1, n-1]
  x[n-1] = b[n-1] / Rt(n-1, n-1);

  // 回代法求解其余方程
  for (int j = n-1; j > 0; --j) {
    double tmp = 0.0;
    for (int k = j; k < n; ++k) {
      tmp += Rt(k, j-1) * x[k];
    }
    x[j-1] = (b[j-1] - tmp) / Rt(j-1, j-1);
  }

  return x;
}