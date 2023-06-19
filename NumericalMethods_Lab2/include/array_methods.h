#pragma once
template<class T>
T** create_array(size_t sizeY, size_t sizeX) {
  T* arr = new T[sizeX * sizeY];
  T** arr2 = new T*[sizeY];
  for (int i = 0; i < sizeY; ++i) {
    arr2[i] = &arr[i * sizeX];
  }
  return arr2;
}

template <class T>
void delete_array(T**& arr) {
  delete (*arr);
  delete arr;
}

template <class T>
void copy_array(T** out, T** in, int sizeY, int sizeX) {
  for (int j = 0; j < sizeY; ++j)
    for (int i = 0; i < sizeX; ++i) in[j][i] = out[j][i];
}

template <class T>
void copy_array(T** out, std::vector<std::vector<T>>& in, int sizeY, int sizeX) {
  in.resize(sizeY);
  for (int j = 0; j < sizeY; ++j) {
    in[j].resize(sizeX);
    for (int i = 0; i < sizeX; ++i) in[j][i] = out[j][i];
  }
}

inline void fill_array_func(double** arr, int** mask, int n, int m, double a, double c,
                     double h, double k, double (*func)(double, double)) {
  double x, y;
  for (int j = 0; j <= m; ++j) {
    y = c + j * k;
    for (int i = 0; i <= n; ++i) {
      if (mask[j][i] == 2 || mask[j][i] == 1) {
        x = a + h * i;
        arr[j][i] = func(x, y);
      }
    }
  }
}

inline int** create_mask_area(int n, int m, int area_x, int area_y, int** area) {
  int** arr2 = create_array<int>(m + 1, n + 1);

  int dx = n / area_x;
  int dy = m / area_y;

  for (int j = 0; j < m + 1; ++j)
    for (int i = 0; i < n + 1; ++i) arr2[j][i] = 0;

  for (int sh_x = 0; sh_x <= area_x; sh_x++) {
    int st_x = sh_x * dx;
    for (int sh_y = 0; sh_y < area_y; sh_y++) {
      int st_y = sh_y * dy;
      bool z, l;
      z = !((sh_x - 1 < 0) || area[sh_y][sh_x - 1] == '0');
      l = !((sh_x >= area_x) || area[sh_y][sh_x] == '0');
      bool flag = z xor l;
      if (flag) {
        for (int b = st_y; b <= st_y + dy; ++b) {
          arr2[b][st_x] = 1;
        }
      }
    }
  }
  for (int sh_y = 0; sh_y <= area_y; sh_y++) {
    int st_y = sh_y * dy;
    for (int sh_x = 0; sh_x < area_x; sh_x++) {
      int st_x = sh_x * dx;
      bool z, l;
      z = !((sh_y - 1 < 0) || area[sh_y - 1][sh_x] == '0');
      l = !((sh_y >= area_y) || area[sh_y][sh_x] == '0');
      bool flag = z xor l;
      if (flag) {
        for (int b = st_x; b <= st_x + dx; ++b) {
          arr2[st_y][b] = 1;
        }
      }
    }
  }
  int tmp;
  for (int y = 0; y <= m; ++y) {
    int flag = 0;
    for (int x = 0; x <= n; ++x) {
      tmp = arr2[y][x];
      switch (tmp) {
        case 1: {
          while (x <= n && arr2[y][x] == 1) ++x;
          if (x <= n) {
            int _x = x / dx;
            int _y = y / dy;
            if (_x >= area_x || _y >= area_y)
              flag = 0;
            else {
              flag = area[_y][_x] == '0' ? 0 : 2;
            }
            arr2[y][x] = flag;
          }
          break;
        }
        case 0: {
          arr2[y][x] = flag;
          break;
        }
      }
    }
  }
  return arr2;
}

inline void fill_border_of_array(double** arr, int** mask, const global_data& data) {
  double a = data.a, b = data.b, c = data.c, d = data.d;
  double h = (data.b - data.a) / data.n;
  double k = (data.d - data.c) / data.m;
  double (*_m)(double, double) = data.func_m;
  double x, y;
  double f1, f2;

  for (int j = 0; j <= data.m; ++j) {
    y = c + j * k;
    f1 = _m(a, y), f2 = _m(b, y);
    for (int i = 0; i <= data.n; ++i) {
      x = a + i * h;
      switch (mask[j][i]) {
        case 0:
          arr[j][i] = 0;
          break;
        case 1:
          arr[j][i] = data.func_m(x, y);
          break;
        case 2:
          // линейна€ интерпл€ци€ значений по оси X
          arr[j][i] = (f2 * (x - a) - f1 * (x - b)) / (b - a);
          break;
      }
    }
  }
}

inline void fill_arrays_of_coord(global_data& data) {
  data.coord_x.resize(data.n + 1);
  data.coord_y.resize(data.m + 1);

  double a = data.a;
  double c = data.c;
  double h = (data.b - data.a) / data.n;
  double k = (data.d - data.c) / data.m;

  for (int j = 0; j <= data.m; ++j) data.coord_y[j] = c + j * k;
  for (int i = 0; i <= data.n; ++i) data.coord_x[i] = a + i * h;
}

inline std::pair<double**, int**> create_area(const global_data& data) {
  std::pair<double**, int**> ans;
  double** arr = create_array<double>(data.m + 1, data.n + 1);
  int** mask =
      create_mask_area(data.n, data.m, data.area_x, data.area_y, data.area);
  fill_border_of_array(arr, mask, data);

  ans.first = arr;
  ans.second = mask;
  return ans;
}

inline void delete_area(std::pair<double**, int**> area) {
  delete_array(area.first);
  delete_array(area.second);
}

template <class T>
void show_array(std::ostream& out, T** arr, int sizeY, int sizeX) {
  for (int j = 0; j < sizeY; ++j) {
    for (int i = 0; i < sizeX; ++i) {
      out << arr[j][i] << " ";
    }
    out << "\n";
  }
}


/* template <class T>
std::pair<int, int> compress_array(T** arr1, int sizeY, int sizeX, T** arr2,
int& new_sizeY, int& new_sizeX, int k) { new_sizeY = sizeY / k; new_sizeX =
sizeX / k; arr2 = create_array<T>(new_sizeY, new_sizeX); for (int j = 0; j <
new_sizeY; ++j){ for (int i = 0; i < new_sizeX; ++i) { arr2[j][i] = arr1[j *
k][i * k];
        }
    }
}*/