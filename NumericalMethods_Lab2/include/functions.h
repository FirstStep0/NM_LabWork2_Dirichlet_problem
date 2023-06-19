#pragma once
#include "../include/structures.h";

double m1(double x, double y) {
  if (abs(x - 1.0) < 1e-6)
    return (y - 2) * (y - 3);
  else if (abs(x - 2.0) < 1e-6)
    return y * (y - 2) * (y - 3);
  else if (abs(y - 2.0) < 1e-6)
    return (x - 1) * (x - 2);
  else if (abs(y - 3.0) < 1e-6)
    return x * (x - 1) * (x - 2);
  else {
    throw std::exception("m1 for other area");
    return 0.0;
  }
}

double m1_test(double x, double y) { return sin(PI * x * y); }

double f1(double x, double y) { return -exp(-x * y * y); }

double f1_test(double x, double y) {
  return PI * PI * sin(PI * x * y) * (x * x + y * y);
}

double u1_test(double x, double y) { return sin(PI * x * y); }

void set_function(int numberVariant, bool test, global_data& data) {
    switch (numberVariant) {
    case 1:
      if (test) {
        data.func_m = m1_test;
        data.func_f = f1_test;
        data.func_u = u1_test;
      } else {
        data.func_m = m1;
        data.func_f = f1;
      }
      break;
    default:
      throw std::exception("incorrect numberVariant");
    }
}