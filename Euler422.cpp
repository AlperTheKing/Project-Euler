#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using namespace std;

namespace {

constexpr int64_t MOD = 1'000'000'007LL;
constexpr int64_t MOD_EXP = MOD - 1;

int64_t mul_mod(int64_t a, int64_t b, int64_t mod) {
  return static_cast<int64_t>((__int128)a * b % mod);
}

int64_t mod_pow(int64_t base, int64_t exp, int64_t mod) {
  int64_t result = 1 % mod;
  int64_t cur = (base % mod + mod) % mod;
  while (exp > 0) {
    if (exp & 1LL) {
      result = mul_mod(result, cur, mod);
    }
    cur = mul_mod(cur, cur, mod);
    exp >>= 1LL;
  }
  return result;
}

int64_t mod_inv(int64_t x) { return mod_pow(x, MOD - 2, MOD); }

pair<int64_t, int64_t> fib_pair_mod(uint64_t n, int64_t mod) {
  if (n == 0) {
    return {0, 1 % mod};
  }
  auto [a, b] = fib_pair_mod(n >> 1, mod); // a=F(k), b=F(k+1)
  int64_t two_b_minus_a = (2 * b - a) % mod;
  if (two_b_minus_a < 0) {
    two_b_minus_a += mod;
  }
  int64_t c = mul_mod(a, two_b_minus_a, mod);                 // F(2k)
  int64_t d = (mul_mod(a, a, mod) + mul_mod(b, b, mod)) % mod; // F(2k+1)
  if (n & 1ULL) {
    return {d, (c + d) % mod};
  }
  return {c, d};
}

uint64_t pow_u64(uint64_t base, uint32_t exp) {
  uint64_t result = 1;
  while (exp > 0) {
    if (exp & 1U) {
      result *= base;
    }
    base *= base;
    exp >>= 1U;
  }
  return result;
}

struct Rational {
  cpp_int num;
  cpp_int den;

  Rational() : num(0), den(1) {}
  Rational(long long n) : num(n), den(1) {}
  Rational(cpp_int n, cpp_int d) : num(std::move(n)), den(std::move(d)) {
    normalize();
  }

  static cpp_int abs_int(cpp_int v) {
    if (v < 0) {
      v = -v;
    }
    return v;
  }

  static cpp_int gcd_int(cpp_int a, cpp_int b) {
    a = abs_int(a);
    b = abs_int(b);
    while (b != 0) {
      cpp_int t = a % b;
      a = b;
      b = t;
    }
    return a;
  }

  void normalize() {
    if (den == 0) {
      throw runtime_error("Zero denominator in Rational.");
    }
    if (den < 0) {
      num = -num;
      den = -den;
    }
    cpp_int g = gcd_int(num, den);
    num /= g;
    den /= g;
  }
};

Rational operator+(const Rational &a, const Rational &b) {
  return Rational(a.num * b.den + b.num * a.den, a.den * b.den);
}

Rational operator-(const Rational &a, const Rational &b) {
  return Rational(a.num * b.den - b.num * a.den, a.den * b.den);
}

Rational operator*(const Rational &a, const Rational &b) {
  return Rational(a.num * b.num, a.den * b.den);
}

Rational operator/(const Rational &a, const Rational &b) {
  return Rational(a.num * b.den, a.den * b.num);
}

struct ExactPoint {
  Rational x;
  Rational y;
};

Rational exact_u_n(uint64_t n) {
  if (n == 1) {
    return Rational(100);
  }
  if (n == 2) {
    return Rational(-75, 2);
  }

  Rational u_im2(100);
  Rational u_im1(cpp_int(-75), cpp_int(2));
  Rational u_i;
  for (uint64_t i = 3; i <= n; ++i) {
    u_i = Rational(25) * u_im2 / u_im1; // u_i = 25*u_{i-2}/u_{i-1}
    u_im2 = u_im1;
    u_im1 = u_i;
  }
  return u_im1;
}

ExactPoint exact_point(uint64_t n) {
  Rational u = exact_u_n(n);
  Rational uu = u * u;
  Rational x = (Rational(3) * uu + Rational(2500)) / (Rational(25) * u);
  Rational y = (Rational(4) * uu - Rational(1875)) / (Rational(25) * u);
  return {x, y};
}

int64_t cpp_mod(const cpp_int &v, int64_t mod) {
  cpp_int r = v % mod;
  if (r < 0) {
    r += mod;
  }
  return static_cast<int64_t>(r);
}

int64_t exact_checksum_mod(uint64_t n) {
  ExactPoint p = exact_point(n);
  int64_t a = cpp_mod(p.x.num, MOD);
  int64_t b = cpp_mod(p.x.den, MOD);
  int64_t c = cpp_mod(p.y.num, MOD);
  int64_t d = cpp_mod(p.y.den, MOD);
  return (a + b + c + d) % MOD;
}

int64_t solve_checksum_mod(uint64_t n) {
  if (n == 0) {
    throw runtime_error("n must be >= 1");
  }

  const uint64_t k = n - 1;
  auto [fk_mod, fk1_mod] = fib_pair_mod(k, MOD_EXP); // F_k, F_{k+1}
  const int64_t F = fk_mod;
  int64_t L = (2 * fk1_mod - fk_mod) % MOD_EXP; // Lucas_k
  if (L < 0) {
    L += MOD_EXP;
  }

  const int sign = (k % 3ULL == 0ULL) ? +1 : -1; // (-1)^{F_k}

  const int64_t two_pow_L = mod_pow(2, L, MOD);
  const int64_t three_pow_F = mod_pow(3, F, MOD);

  int64_t alpha;
  int64_t beta;
  int64_t g1;
  int64_t g2;

  if (n & 1ULL) {
    // odd n: alpha=2^L, beta=3^F
    alpha = two_pow_L;
    beta = three_pow_F;
    g1 = (n == 1 ? 4 : 12);
    g2 = 1;
  } else {
    // even n: alpha=3^F, beta=2^L
    alpha = three_pow_F;
    beta = two_pow_L;
    g1 = 1;
    g2 = (n == 2 ? 6 : 12);
  }

  const int64_t alpha2 = mul_mod(alpha, alpha, MOD);
  const int64_t beta2 = mul_mod(beta, beta, MOD);

  int64_t n1 = (3 * alpha2 + 4 * beta2) % MOD; // for x numerator (pre-reduction)
  int64_t n2 = (4 * alpha2 - 3 * beta2) % MOD; // for y numerator (pre-reduction)
  if (n2 < 0) {
    n2 += MOD;
  }

  const int64_t d_base = mul_mod(two_pow_L, three_pow_F, MOD); // 2^L * 3^F

  int64_t a = mul_mod(n1, mod_inv(g1), MOD);
  int64_t b = mul_mod(d_base, mod_inv(g1), MOD);
  int64_t c = mul_mod(n2, mod_inv(g2), MOD);
  int64_t d = mul_mod(d_base, mod_inv(g2), MOD);

  if (sign < 0) {
    a = (MOD - a) % MOD;
    c = (MOD - c) % MOD;
  }

  return (a + b + c + d) % MOD;
}

bool same_fraction(const Rational &r, long long num, long long den) {
  return r.num == num && r.den == den;
}

bool validate_given_points() {
  struct Check {
    uint64_t n;
    long long x_num;
    long long x_den;
    long long y_num;
    long long y_den;
  };

  const vector<Check> checks = {
      {3, -19, 2, -229, 24},
      {4, 1267, 144, -37, 12},
      {7, 17194218091LL, 143327232LL, 274748766781LL, 1719926784LL},
  };

  for (const auto &ch : checks) {
    ExactPoint p = exact_point(ch.n);
    if (!same_fraction(p.x, ch.x_num, ch.x_den) ||
        !same_fraction(p.y, ch.y_num, ch.y_den)) {
      cerr << "Checkpoint failed at n=" << ch.n << "\n";
      cerr << "Expected x=" << ch.x_num << "/" << ch.x_den << ", y=" << ch.y_num
           << "/" << ch.y_den << "\n";
      cerr << "Got x=" << p.x.num << "/" << p.x.den << ", y=" << p.y.num << "/"
           << p.y.den << "\n";
      return false;
    }
  }

  if (solve_checksum_mod(7) != 806236837LL) {
    cerr << "Checkpoint failed: checksum(n=7) mismatch.\n";
    return false;
  }

  return true;
}

bool validate_parallel_exact_vs_fast(int max_n) {
  atomic<int> next_n{1};
  atomic<bool> ok{true};
  mutex err_mutex;
  vector<string> errors;

  int thread_count = static_cast<int>(thread::hardware_concurrency());
  if (thread_count <= 0) {
    thread_count = 4;
  }
  thread_count = min(thread_count, max_n);
  thread_count = max(thread_count, 1);

  auto worker = [&]() {
    while (true) {
      int n = next_n.fetch_add(1);
      if (n > max_n) {
        break;
      }
      int64_t fast_val = solve_checksum_mod(static_cast<uint64_t>(n));
      int64_t exact_val = exact_checksum_mod(static_cast<uint64_t>(n));
      if (fast_val != exact_val) {
        ok.store(false);
        lock_guard<mutex> lock(err_mutex);
        errors.push_back("Mismatch at n=" + to_string(n) + ": fast=" +
                         to_string(fast_val) + ", exact=" +
                         to_string(exact_val));
      }
    }
  };

  vector<thread> threads;
  threads.reserve(thread_count);
  for (int i = 0; i < thread_count; ++i) {
    threads.emplace_back(worker);
  }
  for (auto &t : threads) {
    t.join();
  }

  if (!ok.load()) {
    for (const string &s : errors) {
      cerr << s << "\n";
    }
    return false;
  }
  return true;
}

bool validate_all() {
  if (!validate_given_points()) {
    return false;
  }
  if (!validate_parallel_exact_vs_fast(15)) {
    return false;
  }
  return true;
}

} // namespace

int main(int argc, char **argv) {
  uint64_t n = pow_u64(11, 14);
  bool run_validation = true;

  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if (arg == "--no-validate") {
      run_validation = false;
    } else {
      n = stoull(arg);
    }
  }

  if (run_validation) {
    if (!validate_all()) {
      return 1;
    }
    cout << "Validation passed (P3, P4, P7, checksum(7), and n<=15 cross-check)."
         << "\n";
  }

  cout << solve_checksum_mod(n) << "\n";
  return 0;
}

