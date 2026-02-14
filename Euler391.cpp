#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

class Solver {
public:
  static constexpr int MAX_N = 1000;

  Solver() : table_a(STRIDE * STRIDE, 0), table_b(STRIDE * STRIDE, 0) {
    a = table_a.data();
    b = table_b.data();
    for (int v = 0; v <= MAX_N; ++v) {
      a[index(0, v)] = v;
      b[index(0, v)] = v;
    }
  }

  int compute_m(int n) {
    const int step = n + 1;
    for (int s = 1; s <= step; ++s) {
      for (int m = 1; m <= s; ++m) {
        const int add = s - m + 1;
        const int prev_row = (m - 1) << SHIFT;
        const int cur_row = m << SHIFT;
        for (int v = 0; v <= n; ++v) {
          int value = b[prev_row + v] + add;
          if (value >= step) {
            value = 0;
          }
          b[cur_row + v] = a[prev_row + value];
        }
      }

      const int row = s << SHIFT;
      const int value = b[row];
      int v = 1;
      while (v <= n && b[row + v] == value) {
        ++v;
      }
      if (v > n) {
        return value;
      }
      swap(a, b);
    }
    return a[index(step, 0)];
  }

private:
  static constexpr int SHIFT = 10;
  static constexpr int STRIDE = 1 << SHIFT;

  vector<int> table_a;
  vector<int> table_b;
  int *a;
  int *b;

  static constexpr int index(int m, int v) { return (m << SHIFT) + v; }
};

int main() {
  Solver solver;
  int m2 = -1;
  int m7 = -1;
  int m20 = -1;
  int64_t sum = 0;

  for (int n = 1; n <= Solver::MAX_N; ++n) {
    const int m = solver.compute_m(n);
    if (n == 2) {
      m2 = m;
    }
    if (n == 7) {
      m7 = m;
    }
    if (n == 20) {
      m20 = m;
    }
    sum += static_cast<int64_t>(m) * m * m;

    if (n == 20) {
      if (m2 != 2 || m7 != 1 || m20 != 4 || sum != 8150) {
        cerr << "Validation failed.\n";
        return 1;
      }
    }
  }

  cout << sum << '\n';
  return 0;
}
