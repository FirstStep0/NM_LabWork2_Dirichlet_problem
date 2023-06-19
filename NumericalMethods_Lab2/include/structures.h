#pragma once
#include<vector>
#include "../include/nlohmann/json.hpp"
using json = nlohmann::json;

struct global_data {
  int numberVariant;
  bool test;

  // stop condition
  std::vector<int> nmax;
  std::vector<double> eps;

  // area
  int n, m;
  double a, b, c, d;
  int** area;
  int area_x, area_y;

  // method
  int numberMethod;
  std::vector<double> param;

  //functions
  double (*func_m)(double, double);
  double (*func_f)(double, double);
  double (*func_u)(double, double);

  //part_of_answer
  double err = LDBL_MIN;
  double x_err = LDBL_MIN, y_err = LDBL_MIN;
  std::vector<std::vector<std::vector<double>>> arr_u;
  std::vector<std::vector<double>> arr_err;
  std::vector<int> N;
  std::vector<double> acc;
  std::vector<double> R;

  // auxiliary
  std::vector<double> coord_x;
  std::vector<double> coord_y;

  void initialize() {
    if (test) {
      N.resize(1);
      acc.resize(1);
      R.resize(1);
    } else {
      N.resize(2);
      acc.resize(2);
      R.resize(2);
    }
    arr_u.resize(2);
    arr_err.resize(1);

    double h, k;
    h = (b - a) / n;
    k = (d - c) / m;

    //for (int i = 0; i <= n; ++i) coord_x[i] = a + i * h;
    //for (int j = 0; j <= m; ++j) coord_y[j] = c + j * k;
  }

  inline std::string get_answer() {
    json j;
    j["a"] = a;
    j["b"] = b;
    j["c"] = c;
    j["d"] = d;
    j["n"] = n;
    j["m"] = m;
    j["N"] = N;
    j["nmax"] = nmax;
    j["eps"] = eps;
    j["test"] = test;
    j["err"] = err;
    j["acc"] = acc;
    j["R"] = R;
    j["x_err"] = x_err;
    j["y_err"] = y_err;
    j["param"] = param;
    j["arr_u"] = arr_u;
    j["arr_err"] = arr_err;
    return j.dump();
  };
};

struct result_method {
  int count;
  double R;
  double acc;
};