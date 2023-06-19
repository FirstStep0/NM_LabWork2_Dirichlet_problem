#include "../include/ConjugateGradientsMethod.h"
#include "../include/methods.h"
#include <valarray>
#include <omp.h>

static double H, K, A;
static int _n, _m;

static std::valarray<double> mult_by_A(const std::valarray<double>& x)
{
    std::valarray<double> res = x;

    #pragma omp parallel for
    for(int i = 1; i < _n - 1; i++)
        for (int j = 1; j < _m - 1; j++)
        {
            double tmp = 0;
            tmp = A * x[i * _m + j];
            if (i - 1 > 0)
                tmp += K * x[(i - 1) * _m + j];
            if (i + 1 < _n - 1)
                tmp += K * x[(i + 1) * _m + j];

            if (j - 1 > 0)
                tmp += H * x[i * _m + j - 1];
            if (j + 1 < _m - 1)
                tmp += H * x[i * _m + j + 1];

            res[i * _m + j] = tmp;
        }

    return res;
}

static std::valarray<double> mult_by_A_minus_f(const std::valarray<double>& x, const std::valarray<double>& f)
{
    std::valarray<double> res = x;


    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int i = 0; i < _n; i++)
        {
            res[i * _m] = 0;
            res[i * _m + _m - 1] = 0;
        }

        #pragma omp for nowait
        for (int j = 0; j < _m; j++)
        {
            res[j] = 0;
            res[(_n - 1) * _m + j] = 0;
        }

        #pragma omp for nowait
        for (int i = 1; i < _n - 1; i++)
            for (int j = 1; j < _m - 1; j++)
            {
                double tmp = 0;
                tmp = A * x[i * _m + j];
                tmp += K * x[(i - 1) * _m + j];
                tmp += K * x[(i + 1) * _m + j];
                tmp += H * x[i * _m + j - 1];
                tmp += H * x[i * _m + j + 1];
                tmp += f[i * _m + j];

                res[i * _m + j] = tmp;
            }
    }

    return res;
}

static double dot_prod(const std::valarray<double>& x, const std::valarray<double>& y)
{
    double res = 0;
    #pragma omp parallel for reduction(+: res)
    for(int i = 1; i < _n - 1 ; i++)
        for (int j = 1; j < _m - 1; j++)
        {
            res += x[i * _m + j] * y[i * _m + j];
        }

    return res;
}

omp_lock_t mutex;

static double get_norm(const std::valarray<double>& vec)
{
    assert(vec.size());

    double res = 0;

    #pragma omp parallel
    {
        double tmp_res = 0;
        #pragma omp for
        for (int i = 1; i < _n - 1; i++)
            for (int j = 1; j < _m - 1; j++)
            {
                tmp_res = std::max(std::abs(vec[i * _m + j]), tmp_res);
            }

        omp_set_lock(&mutex);
        res = std::max(tmp_res, res);
        omp_unset_lock(&mutex);
    }

    return res;
}

result_method ConjugateGradientsMethod_impl(double** v, double** f, int** mask,
    int n,
    int m, double h, double k, int nmax,
    double eps,
    const std::vector<double>& param) {

    _n = m + 1;
    _m = n + 1;

    H = 1.0 / (h * h);
    K = 1.0 / (k * k);
    A = -2.0 * (H + K);

    std::valarray<double> x(_n * _m);
    std::valarray<double> b(_n * _m);
    std::copy(*v, *v + _n * _m, std::begin(x));
    std::copy(*f, *f + _n * _m, std::begin(b));

    omp_init_lock(&mutex);

    // First step
    std::valarray<double> h0;
    std::valarray<double> A_h;
    double dot_prod_A_h_h;
    double alpha;

    // Other step
    double acc = std::numeric_limits<double>::max();
    uint64_t count = 0;
    // nmax = std::min(nmax, n * m);
    while (count < nmax && acc > eps) {

        if (count++ % 10000000 == 0)
        {
            h0 = -1 * mult_by_A_minus_f(x, b);
            A_h = mult_by_A(h0);
            dot_prod_A_h_h = dot_prod(A_h, h0);
            alpha = -1 * dot_prod(mult_by_A_minus_f(x, b), h0) / dot_prod_A_h_h;
            x += alpha * h0;
        }
        else
        {
            auto r1 = mult_by_A_minus_f(x, b);
            double beta = dot_prod(A_h, r1) / dot_prod_A_h_h;
            h0 = -1 * r1 + beta * h0;
            A_h = mult_by_A(h0);
            dot_prod_A_h_h = dot_prod(A_h, h0);
            alpha = -1 * dot_prod(r1, h0) / dot_prod_A_h_h;
            x += alpha * h0;
        }
        acc = get_norm(alpha * h0);
    }

    std::copy(std::begin(x), std::end(x), *v);

    result_method res;
    res.count = count;
    res.R = get_norm(mult_by_A_minus_f(x, b));
    res.acc = acc;

    omp_destroy_lock(&mutex);

    return res;
}