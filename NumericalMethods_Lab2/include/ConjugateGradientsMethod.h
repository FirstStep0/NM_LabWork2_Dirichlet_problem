#pragma once
#include <vector>
#include "../include/structures.h"


result_method ConjugateGradientsMethod_impl(double** v, double** f, int** mask,
    int n,
    int m, double h, double k, int nmax,
    double eps,
    const std::vector<double>& param);