#pragma once
#include "../include/structures.h"
#include "../include/array_methods.h"
#include "../include/ConjugateGradientsMethod.h"

/* void updateErr(double err, double x, double y, answer& ans) {
  if (err > ans.err) {
    ans.x = x;
    ans.y = y;
    ans.err = err;
  }
}*/
inline double get_R(double** u, double** f, int** mask, int n, int m, double H,
             double K, double A) {
  double R = LDBL_MIN;
  for (int j = 0; j <= m; ++j) {
    for (int i = 0; i <= n; ++i) {
      if (mask[j][i] == 2) {
        R = std::max(R, abs(H * (u[j][i + 1] + u[j][i - 1]) +
                            K * (u[j + 1][i] + u[j - 1][i]) + A * u[j][i] +
                            f[j][i]));  //?
      }
    }
  }
  return R;
}

inline result_method method_upper_relaxation(double** v, double** f, int** mask,
                                      int n, int m,
                             double h, double k, int nmax, double eps,
                             const std::vector<double>& param)
{
  //double** _u = create_array<double>(data.m + 1, data.n + 1);
  //copy_array<double>(u, _u, data.m + 1, data.n + 1);

  double w = param[0];

  double H = 1.0 / (h * h);
  double K = 1.0 / (k * k);
  double A = -2.0 * (H + K);

  double R = LDBL_MIN;
  double acc = LDBL_MAX;
  double x, y, tmp;
  int count = 0;

  while (count < nmax && acc > eps) {
    acc = LDBL_MIN;
    R = LDBL_MIN;
    for (int j = 1; j < m; ++j) {
      for (int i = 1; i < n; ++i) {
        if (mask[j][i] == 2) {
          tmp = -(w * (H * (v[j][i + 1] + v[j][i - 1]) +
                       K * (v[j + 1][i] + v[j - 1][i]) + f[j][i]) +
                  (w - 1.0) * A * v[j][i]) /
                A;
          acc = std::max(acc, abs(v[j][i] - tmp));
          v[j][i] = tmp;
        }
      }
    }
    ++count;
  }
  result_method res;
  res.count = count;
  res.R = get_R(v, f, mask, n, m, H, K, A);
  res.acc = acc;
  return res;
}

inline result_method SimpleIterationMethod(double** v, double** f, int** mask,
                                       int n, int m, double h, double k,
                                       int nmax, double eps,
                                       const std::vector<double>& param) {
  double w = param[0];

  double H = 1.0 / (h * h);
  double K = 1.0 / (k * k);
  double A = -2.0 * (H + K);

  double R = LDBL_MIN;
  double acc = LDBL_MAX;
  double x, y, tmp;
  int count = 0;

  while (count < nmax && acc > eps) {
    acc = LDBL_MIN;
    R = LDBL_MIN;
    for (int j = 1; j < m; ++j) {
      for (int i = 1; i < n; ++i) {
        if (mask[j][i] == 2) {
          //code method
        }
      }
    }
    ++count;
  }
  result_method res;
  res.count = count;
  res.R = get_R(v, f, mask, n, m, H, K, A);
  res.acc = acc;
  return res;
}



inline result_method ConjugateGradientsMethod(double** v, double** f, int** mask,
                                       int n,
                                      int m, double h, double k, int nmax,
                                      double eps,
                                      const std::vector<double>& param) {
    return ConjugateGradientsMethod_impl(v, f, mask, n, m, h, k, nmax, eps, param);
}

inline result_method (*choose_method(int numberMethod))(double** v, double** f, int** mask, int n,
                                 int m, double h, double k, int nmax,
                                 double eps, const std::vector<double>& param) {
  switch (numberMethod){ 
  case 1:
      return method_upper_relaxation;
  case 2:
      return SimpleIterationMethod;
  case 3:
      return ConjugateGradientsMethod;
  default:
      throw std::exception("Method don't exist");
  }
}