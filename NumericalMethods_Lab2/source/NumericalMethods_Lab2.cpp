#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "../include/functions.h"
#include "../include/methods.h"


#define CHECK(x, y) if (!(x)) throw std::exception(y)
#define CHECK_ARRAY(arr, cond, exc) for (int i = 0; i < arr.size(); ++i) { auto item = arr[i];if (!(cond)) throw std::exception(exc);}                                      
int precision = 6;

void parse(std::ifstream& f, global_data& data) {
  json j = json::parse(f);
  data.a = j["a"];
  data.b = j["b"];
  data.c = j["c"];
  data.d = j["d"];
  data.n = j["n"];
  data.m = j["m"];
  data.nmax = j["nmax"].get<std::vector<int>>();
  data.eps = j["eps"].get<std::vector<double>>();
  data.test = j["test"];
  data.numberMethod = j["numberMethod"];
  data.numberVariant = j["numberVariant"];
  data.param = j["param"].get<std::vector<double>>();
  data.area_x = j["area_x"];
  data.area_y = j["area_y"];
  std::string str = j["area"];

  CHECK(data.a < data.b, "a >= b");
  CHECK(data.c < data.d, "c >= d");
  CHECK(data.n >= 0, "n < 0");
  CHECK(data.m >= 0, "m < 0");
  CHECK(data.numberMethod >= 0, "numberMethod < 0");
  CHECK(data.numberVariant >= 0, "numberVariant < 0");
  CHECK_ARRAY(data.nmax, item >= 0, "nmax < 0");
  CHECK_ARRAY(data.eps, item >= 0, "eps < 0");
  CHECK(data.area_x > 0, "area_x <= 0");
  CHECK(data.area_y > 0, "area_y <= 0");
  CHECK(data.n % data.area_x == 0, "n % area_x != 0");
  CHECK(data.m % data.area_y == 0, "m % area_y != 0");

  data.area = create_array<int>(data.area_y, data.area_x);
  for (int j = 0; j < data.area_y; ++j) {
    for (int i = 0; i < data.area_x; ++i) {
      data.area[j][i] = str[j * data.area_x + i];
    }
  }

  data.initialize();
}

void add_result_method(int pos, double** v, result_method& res, global_data& data) {
    if (data.arr_u.size() <= pos) data.arr_u.resize(pos + 1);
    if (data.acc.size() <= pos) data.acc.resize(pos + 1);
    if (data.R.size() <= pos) data.R.resize(pos + 1);
    if (data.N.size() <= pos) data.N.resize(pos + 1);
    copy_array(v, data.arr_u[pos], data.m + 1, data.n + 1);
    data.acc[pos] = res.acc;
    data.R[pos] = res.R;
    data.N[pos] = res.count;
}

void set_error(std::vector<std::vector<double>>& arr1,
               std::vector<std::vector<double>>& arr2,
               std::vector<std::vector<double>>& arr_out, int** mask,
               global_data& data, int step) {
    double& err = data.err;
    err = LDBL_MIN;
    int _i, _j;
    _i = _j = -1;
    double tmp;
    arr_out.resize(data.m + 1);
    for (int j = 0; j <= data.m; ++j) {
        arr_out[j].resize(data.n + 1);
        for (int i = 0; i <= data.n; ++i) {
            if (mask[j][i] == 2 || mask[j][i] == 1) {
            tmp = abs(arr2[j * step][i * step] - arr1[j][i]);
            arr_out[j][i] = tmp;
            if (err < tmp) {
                err = tmp;
                _i = i;
                _j = j;
            }
            }
        }
    }
    data.x_err = data.a + _i * (data.b - data.a) / data.n;
    data.y_err = data.c + _j * (data.d - data.c) / data.m;
}

int main(int argc, char* argv[]) {
    std::cout.precision(precision);

    clock_t start = clock();

    std::ifstream in;
    std::ofstream out;
    global_data data;

    try {
        CHECK(argc == 3, "input argument as [input path] [output path]");

        in.open(argv[1]);
        if (!in.is_open()) throw std::exception("path not open");
        parse(in, data);
        in.close();
        set_function(data.numberVariant, data.test, data);

        double h, k;
        h = (data.b - data.a) / data.n;
        k = (data.d - data.c) / data.m;

        auto p = create_area(data);
        double** v = p.first;
        int** mask = p.second;
        // show_array(std::cout, v, data.m + 1, data.n + 1);
        double** f = create_array<double>(data.m + 1, data.n + 1);
        fill_array_func(f, mask, data.n, data.m, data.a, data.c, h, k,
                        data.func_f);

        auto method = choose_method(data.numberMethod);

        auto& coord_x = data.coord_x;
        auto& coord_y = data.coord_y;
        double x, y;

        if (data.test) {
            result_method res =
                method(v, f, mask, data.n, data.m, h, k,
                                        data.nmax[0], data.eps[0], data.param);
            double** u = create_array<double>(data.m + 1, data.n + 1);

            fill_array_func(u, mask, data.n, data.m, data.a, data.c, h, k,
                            data.func_u);
            add_result_method(0, v, res, data);
            copy_array(u, data.arr_u[1], data.m + 1, data.n + 1);
            set_error(data.arr_u[0], data.arr_u[1], data.arr_err, mask, data,
                      1);
        } else {
            result_method res =
                method(v, f, mask, data.n, data.m, h, k,
                                        data.nmax[0], data.eps[0], data.param);
            add_result_method(0, v, res, data);

            // task 2
            data.n *= 2;
            data.m *= 2;

            auto p2 = create_area(data);
            double** v2 = p2.first;
            int** mask2 = p2.second;
            double** f2 = create_array<double>(data.m + 1, data.n + 1);

            h = (data.b - data.a) / data.n;
            k = (data.d - data.c) / data.m;

            fill_array_func(f2, mask2, data.n, data.m, data.a, data.c, h, k,
                            data.func_f);
            res =
                method(v2, f2, mask2, data.n, data.m, h, k,
                                        data.nmax[1], data.eps[1], data.param);
            add_result_method(1, v2, res, data);
            data.m /= 2;
            data.n /= 2;
            set_error(data.arr_u[0], data.arr_u[1], data.arr_err, mask, data,
                      2);
        }
        out.open(argv[2]);
        out.clear();
        out << data.get_answer() << "\n";
        out.close();

        delete_area(p);
        delete_array(f);
    }

    catch (const std::exception& e) {
        std::cerr << e.what();
        return -1;
    }

    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    std::cout << "The time: " << seconds << " seconds\n";
}
     