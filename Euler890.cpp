#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

constexpr int64_t kMod = 1'000'000'007;

std::vector<int> ToBitsLSB(cpp_int n) {
  std::vector<int> bits;
  if (n == 0) {
    return bits;
  }
  while (n > 0) {
    bits.push_back(static_cast<int>(n & 1));
    n >>= 1;
  }
  return bits;
}

int64_t ModAdd(int64_t a, int64_t b) {
  int64_t s = a + b;
  if (s >= kMod) {
    s -= kMod;
  }
  return s;
}

int64_t ModMul(int64_t a, int64_t b) {
  return static_cast<int64_t>((__int128)a * b % kMod);
}

int64_t ComputePartitions(const cpp_int &n, int thread_count = 1) {
  if (n == 0) {
    return 1;
  }

  std::vector<int> bits = ToBitsLSB(n);
  if (bits.size() == 1) {
    return 1;
  }

  const int L = static_cast<int>(bits.size());

  std::vector<int64_t> pow2(L + 3, 1);
  for (int i = 1; i < static_cast<int>(pow2.size()); ++i) {
    pow2[i] = ModAdd(pow2[i - 1], pow2[i - 1]);
  }

  std::vector<std::vector<int64_t>> C(L + 1, std::vector<int64_t>(L + 1, 0));
  for (int i = 0; i <= L; ++i) {
    C[i][0] = 1;
    for (int j = 1; j <= i; ++j) {
      C[i][j] = ModAdd(C[i - 1][j - 1], C[i - 1][j]);
    }
  }

  std::vector<std::vector<int64_t>> w(L + 1);
  std::vector<std::vector<int64_t>> w2(L + 1);
  for (int j = 0; j <= L; ++j) {
    w[j].resize(j + 1);
    for (int m = 0; m <= j; ++m) {
      w[j][m] = ModMul(C[j][m], pow2[j - m]);
    }
    w2[j].resize(j + 2);
    w2[j][0] = w[j][0];
    for (int m = 1; m <= j; ++m) {
      int64_t val = ModAdd(w[j][m], w[j][m - 1]);
      w2[j][m] = val;
    }
    w2[j][j + 1] = w[j][j];
  }

  std::vector<std::vector<int64_t>>().swap(C);

  std::vector<int64_t> a(1, 1);

  auto compute_range = [&](const std::vector<int64_t> &b,
                           std::vector<int64_t> &next, int bit, int d, int start,
                           int end) {
    for (int j = start; j < end; ++j) {
      int max_m = 0;
      if (bit == 0) {
        max_m = std::min(j, d + 1 - j);
      } else {
        max_m = std::min(j + 1, d + 1 - j);
      }
      int64_t sum = 0;
      if (bit == 0) {
        const auto &wj = w[j];
        for (int m = 0; m <= max_m; ++m) {
          int i = j + m;
          sum += ModMul(b[i], wj[m]);
          if (sum >= kMod) {
            sum -= kMod;
          }
        }
      } else {
        const auto &w2j = w2[j];
        for (int m = 0; m <= max_m; ++m) {
          int i = j + m;
          sum += ModMul(b[i], w2j[m]);
          if (sum >= kMod) {
            sum -= kMod;
          }
        }
      }
      next[j] = sum;
    }
  };

  const int min_parallel = 128;
  if (thread_count < 1) {
    thread_count = 1;
  }

  for (int idx = 1; idx < L; ++idx) {
    const int bit = bits[idx];
    const int d = static_cast<int>(a.size()) - 1;

    std::vector<int64_t> b(d + 2, 0);
    b[0] = a[0];
    for (int i = 1; i <= d; ++i) {
      b[i] = ModAdd(a[i], a[i - 1]);
    }
    b[d + 1] = a[d];

    std::vector<int64_t> next(d + 2, 0);

    if (thread_count == 1 || d + 2 < min_parallel) {
      compute_range(b, next, bit, d, 0, d + 2);
    } else {
      const int threads = std::min(thread_count, d + 2);
      const int chunk = (d + 2 + threads - 1) / threads;
      std::vector<std::thread> pool;
      pool.reserve(threads);
      for (int t = 0; t < threads; ++t) {
        int start = t * chunk;
        int end = std::min(d + 2, start + chunk);
        if (start >= end) {
          continue;
        }
        pool.emplace_back(compute_range, std::cref(b), std::ref(next), bit, d,
                          start, end);
      }
      for (auto &th : pool) {
        th.join();
      }
    }

    a.swap(next);
  }

  return a[0] % kMod;
}

bool Validate() {
  cpp_int n7 = 7;
  if (ComputePartitions(n7) != 6) {
    std::cerr << "Validation failed: p(7) != 6\n";
    return false;
  }

  cpp_int n7_7 = 1;
  for (int i = 0; i < 7; ++i) {
    n7_7 *= 7;
  }
  if (ComputePartitions(n7_7) != 144548435) {
    std::cerr << "Validation failed: p(7^7) != 144548435\n";
    return false;
  }

  return true;
}

}  // namespace

int main() {
  if (!Validate()) {
    return 1;
  }

  cpp_int n = 1;
  for (int i = 0; i < 777; ++i) {
    n *= 7;
  }

  int threads = static_cast<int>(std::thread::hardware_concurrency());
  if (threads <= 0) {
    threads = 1;
  }

  int64_t result = ComputePartitions(n, threads);
  std::cout << result << "\n";
  return 0;
}
