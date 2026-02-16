#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u128 = unsigned __int128;

struct Options {
  bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "--skip-checkpoints") {
      options.run_checkpoints = false;
      continue;
    }
    std::cerr << "Unknown argument: " << arg << '\n';
    return false;
  }
  return true;
}

struct BigInt {
  static constexpr int WORDS = 6;
  u64 data[WORDS];

  BigInt() { std::memset(data, 0, sizeof(data)); }

  explicit BigInt(u64 val) {
    std::memset(data, 0, sizeof(data));
    data[0] = val;
  }

  bool operator<(const BigInt& other) const {
    for (int i = WORDS - 1; i >= 0; --i) {
      if (data[i] != other.data[i]) return data[i] < other.data[i];
    }
    return false;
  }

  bool operator<=(const BigInt& other) const { return !(other < *this); }

  bool operator>(const BigInt& other) const { return other < *this; }

  BigInt operator*(const BigInt& other) const {
    BigInt res;
    for (int i = 0; i < WORDS; ++i) {
      u128 carry = 0;
      for (int j = 0; j < WORDS - i; ++j) {
        const u128 prod = static_cast<u128>(data[i]) * other.data[j] + res.data[i + j] + carry;
        res.data[i + j] = static_cast<u64>(prod);
        carry = prod >> 64;
      }
    }
    return res;
  }

  void multiply_by(u64 p) {
    u128 carry = 0;
    for (int i = 0; i < WORDS; ++i) {
      const u128 prod = static_cast<u128>(data[i]) * p + carry;
      data[i] = static_cast<u64>(prod);
      carry = prod >> 64;
    }
  }

  u64 div_mod(u64 divisor) {
    u128 rem = 0;
    for (int i = WORDS - 1; i >= 0; --i) {
      const u128 cur = data[i] + (rem << 64);
      data[i] = static_cast<u64>(cur / divisor);
      rem = cur % divisor;
    }
    return static_cast<u64>(rem);
  }

  u64 mod_10_16() const {
    BigInt tmp = *this;
    return tmp.div_mod(10000000000000000ULL);
  }

  u128 low_u128() const {
    return static_cast<u128>(data[1]) << 64 | static_cast<u128>(data[0]);
  }
};

bool square_leq(const BigInt& x, const BigInt& bound) {
  const BigInt sq = x * x;
  return sq <= bound;
}

std::vector<u32> primes_below(int limit) {
  std::vector<unsigned char> is_prime(static_cast<std::size_t>(limit), 1);
  is_prime[0] = 0;
  is_prime[1] = 0;
  for (int p = 2; p * p < limit; ++p) {
    if (!is_prime[static_cast<std::size_t>(p)]) continue;
    for (int q = p * p; q < limit; q += p) is_prime[static_cast<std::size_t>(q)] = 0;
  }
  std::vector<u32> primes;
  for (int p = 2; p < limit; ++p) {
    if (is_prime[static_cast<std::size_t>(p)]) primes.push_back(static_cast<u32>(p));
  }
  return primes;
}

struct Subset {
  double log_value;
  u32 mask;

  bool operator<(const Subset& other) const { return log_value < other.log_value; }
};

std::vector<Subset> enumerate_subsets(const std::vector<u32>& primes) {
  const std::size_t n = primes.size();
  const std::size_t m = std::size_t{1} << n;
  std::vector<double> logs(n);
  for (std::size_t i = 0; i < n; ++i) logs[i] = std::log(static_cast<double>(primes[i]));

  std::vector<Subset> out(m);
  out[0] = {0.0, 0U};
  for (u32 mask = 1; mask < m; ++mask) {
    const u32 lsb = mask & (~mask + 1U);
    const int bit = __builtin_ctz(lsb);
    const u32 prev = mask ^ lsb;
    out[mask].log_value = out[prev].log_value + logs[static_cast<std::size_t>(bit)];
    out[mask].mask = mask;
  }

  std::sort(out.begin(), out.end());
  return out;
}

BigInt build_from_masks(const std::vector<u32>& left_primes, const std::vector<u32>& right_primes,
                        u32 left_mask, u32 right_mask) {
  BigInt v(1);
  u32 lm = left_mask;
  while (lm) {
    const int b = __builtin_ctz(lm);
    v.multiply_by(static_cast<u64>(left_primes[static_cast<std::size_t>(b)]));
    lm &= (lm - 1U);
  }
  u32 rm = right_mask;
  while (rm) {
    const int b = __builtin_ctz(rm);
    v.multiply_by(static_cast<u64>(right_primes[static_cast<std::size_t>(b)]));
    rm &= (rm - 1U);
  }
  return v;
}

BigInt solve_for_primes(const std::vector<u32>& primes) {
  BigInt product_all(1);
  double total_log = 0.0;
  for (u32 p : primes) {
    product_all.multiply_by(static_cast<u64>(p));
    total_log += std::log(static_cast<double>(p));
  }
  const double target_log = total_log * 0.5;

  const std::size_t split = primes.size() / 2;
  std::vector<u32> left_primes(primes.begin(), primes.begin() + static_cast<std::ptrdiff_t>(split));
  std::vector<u32> right_primes(primes.begin() + static_cast<std::ptrdiff_t>(split), primes.end());

  const std::vector<Subset> left = enumerate_subsets(left_primes);
  const std::vector<Subset> right = enumerate_subsets(right_primes);

  double best_log = -1.0;
  std::int64_t j = static_cast<std::int64_t>(right.size()) - 1;
  for (const Subset& l : left) {
    while (j > 0 && l.log_value + right[static_cast<std::size_t>(j)].log_value > target_log + 1e-12) --j;
    const double cand = l.log_value + right[static_cast<std::size_t>(j)].log_value;
    if (cand <= target_log + 1e-12 && cand > best_log) best_log = cand;
  }

  BigInt best_value(0);
  j = static_cast<std::int64_t>(right.size()) - 1;
  for (const Subset& l : left) {
    while (j > 0 && l.log_value + right[static_cast<std::size_t>(j)].log_value > target_log + 1e-10) --j;

    for (int t = -1; t <= 12; ++t) {
      const std::int64_t k = j - t;
      if (k < 0 || k >= static_cast<std::int64_t>(right.size())) continue;

      const Subset& r = right[static_cast<std::size_t>(k)];
      const double cand_log = l.log_value + r.log_value;
      if (t >= 0 && cand_log + 1e-9 < best_log) break;
      if (cand_log > target_log + 1e-8) continue;

      const BigInt cand = build_from_masks(left_primes, right_primes, l.mask, r.mask);
      if (square_leq(cand, product_all) && cand > best_value) best_value = cand;
    }
  }

  return best_value;
}

u128 brute_small(const std::vector<u32>& primes) {
  const int n = static_cast<int>(primes.size());
  const u32 all = (n == 32 ? 0xFFFFFFFFU : ((1U << n) - 1U));
  (void)all;

  u128 P = 1;
  for (u32 p : primes) P *= static_cast<u128>(p);

  const u64 lim = (n >= 31) ? 0 : (1ULL << n);
  u128 best = 0;
  for (u64 mask = 0; mask < lim; ++mask) {
    u128 v = 1;
    for (int i = 0; i < n; ++i) {
      if (mask & (1ULL << i)) v *= static_cast<u128>(primes[static_cast<std::size_t>(i)]);
    }
    if (v * v <= P && v > best) best = v;
  }
  return best;
}

bool run_checkpoints() {
  {
    std::vector<u32> p = primes_below(30);
    const BigInt got = solve_for_primes(p);
    if (got.low_u128() != brute_small(p)) {
      std::cerr << "Checkpoint failed for primes<30" << '\n';
      return false;
    }
  }
  {
    std::vector<u32> p = primes_below(40);
    const BigInt got = solve_for_primes(p);
    if (got.low_u128() != brute_small(p)) {
      std::cerr << "Checkpoint failed for primes<40" << '\n';
      return false;
    }
  }
  {
    std::vector<u32> p = primes_below(50);
    const BigInt got = solve_for_primes(p);
    if (got.low_u128() != brute_small(p)) {
      std::cerr << "Checkpoint failed for primes<50" << '\n';
      return false;
    }
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  Options options;
  if (!parse_arguments(argc, argv, options)) return 1;
  if (options.run_checkpoints && !run_checkpoints()) return 2;

  const std::vector<u32> primes = primes_below(190);
  const BigInt ans = solve_for_primes(primes);
  std::cout << ans.mod_10_16() << '\n';
  return 0;
}
