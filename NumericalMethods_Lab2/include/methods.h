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
  double R = -LDBL_MAX;
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

  const double H = 1.0 / (h * h);
  const double K = 1.0 / (k * k);
  const double A = -2.0 * (H + K);

  const double z = -1.0 / A;
  const double l = (w - 1.0) * A;

  double R = -LDBL_MAX;
  double acc = LDBL_MAX;
  double x, y, tmp;
  int count = 0;

  //int left = -sizeof(double);
  const int right = 1;
  const int up = (n + 1);
  int address;

  double* arr = v[0];
  const double* arr_f = f[0];

  while (count < nmax && acc > eps) {
    address = ((n + 1) + 1);
    acc = -LDBL_MAX;
    R = LDBL_MIN;
    for (int j = 1; j < m; ++j) {
      //int address = (j * (n + 1) + 1);
      //std::cout << address << "\n";
      for (int i = 1; i < n; ++i) {
       //if (mask[j][i] == 2) {
          //tmp = -(w * (H * (v[j][i + 1] + v[j][i - 1]) +
                //       K * (v[j + 1][i] + v[j - 1][i]) + f[j][i]) +
              //    (w - 1.0) * A * v[j][i]) /
            //    A;
        tmp = (w * (H * (arr[address - right] + arr[address + right]) +
                    K * (arr[address - up] + arr[address + up]) + arr_f[address]) +
               l * arr[address]) *
              z;
          acc = std::max(acc, abs(arr[address] - tmp));
        arr[address] = tmp;
        //}   
        address += 1;
      }
      address += 2;
    }
    ++count;
    if ((count * 100 - 1) / nmax != (count * 100 / nmax)) {
      std::cout << ((count * 100) / nmax) << "\n";
      std::cout.flush();
    }
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

  // eigenvalue estimates
  double Mmin = 1000000;  // just a big number
  double Mmax = 0;

  for (int j = 1; j < m; ++j) {
      for (int i = 1; i < n; ++i) {
          if (mask[j][i] == 2) {
              double center = A;
              double radius = 0;
              if (mask[j + 1][i] == 2)
                  radius += K;
              if (mask[j - 1][i] == 2)
                  radius += K;
              if (mask[j][i + 1] == 2)
                  radius += H;
              if (mask[j][i - 1] == 2)
                  radius += H;

              if (center - radius < Mmin)
                  Mmin = center - radius;
              if (center + radius > Mmax)
                  Mmax = center + radius;
          }
      }
  }

  double tau = 2 / (Mmin + Mmax);

  // zero approximation
  for (int j = 1; j < m; ++j) {
    for (int i = 1; i < n; ++i) {
      if (mask[j][i] == 2) {
        double b = -f[j][i];
        if (mask[j][i - 1] == 1)
          b -= v[j][i - 1] / H;
        if (mask[j][i + 1] == 1)
          b -= v[j][i + 1] / H;
        if (mask[j - 1][i] == 1)
          b -= v[j - 1][i] / K;
        if (mask[j + 1][i] == 1)
          b -= v[j + 1][i] / K;
        v[j][i] = b / A;
      }
    }
  }

  // simple iteration method
  while (count < nmax && acc > eps) {
    acc = -LDBL_MAX;
    R = LDBL_MIN;
    for (int j = 1; j < m; ++j) {
      for (int i = 1; i < n; ++i) {
        if (mask[j][i] == 2) {
            double b = -f[j][i];
            if (mask[j][i - 1] == 1)
                b += v[j][i - 1] / H;
            if (mask[j][i + 1] == 1)
                b += v[j][i + 1] / H;
            if (mask[j - 1][i] == 1)
                b += v[j - 1][i] / K;
            if (mask[j + 1][i] == 1)
                b += v[j + 1][i] / K;

            double Ax = 0;
            Ax += A * v[j][i];
            if (mask[j + 1][i] == 2)
                Ax += K * v[j + 1][i];
            if (mask[j - 1][i] == 2)
                Ax += K * v[j - 1][i];
            if (mask[j][i + 1] == 2)
                Ax += H * v[j][i + 1];
            if (mask[j][i - 1] == 2)
                Ax += H * v[j][i - 1];

            double R = Ax - b;
            acc = std::max(v[j][i], v[j][i] - tau * R);
            v[j][i] = v[j][i] - tau * R;
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