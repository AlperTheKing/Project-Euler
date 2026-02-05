#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

struct BitsetHash {
  size_t operator()(const vector<uint64_t> &v) const noexcept {
    uint64_t h = 1469598103934665603ULL;
    for (uint64_t w : v) {
      h ^= w;
      h *= 1099511628211ULL;
    }
    return static_cast<size_t>(h);
  }
};

class Solver {
public:
  explicit Solver(int n_val) : n(n_val) {
    B = 1;
    b = 0;
    while (B <= n) {
      B <<= 1;
      ++b;
    }
    level = n + 1 - b;
    wordCount = (B + 63) / 64;

    precompute_popcount();
    precompute_s_base();
    precompute_tmax_in();
    precompute_max_t_out();

    base_cost_to_end.resize(B);
    for (int i = 0; i < B; ++i) {
      base_cost_to_end[i] = s_base[B] - s_base[i + 1];
    }

    pool.emplace_back(wordCount, 0ULL);
    pool_map.emplace(pool[0], 0);

    bit.assign(B + 1, 0);
    status.assign(B, 0);
    prefixLose.assign(B, 0);
  }

  long long compute_M() {
    int id = solve_partial(level, 0, 1, 0);
    const auto &bits = pool[id];

    int s = 0;
    int M = 0;
    for (int k = 1; k < B; ++k) {
      s += popcnt[k];
      if (s > n) {
        break;
      }
      if (!get_bit(bits, k)) {
        M = s;
      }
    }
    return M;
  }

private:
  int n;
  int B;
  int b;
  int level;
  int wordCount;

  vector<int> popcnt;
  vector<int> s_base;
  vector<int> base_cost_to_end;

  vector<vector<uint16_t>> tmax_in;

  vector<vector<int16_t>> max_t_out;

  vector<vector<uint64_t>> pool;
  unordered_map<vector<uint64_t>, int, BitsetHash> pool_map;

  unordered_map<uint64_t, int> memo_full;
  unordered_map<uint64_t, int> memo_partial;

  vector<int> bit;
  vector<uint8_t> status;
  vector<uint8_t> prefixLose;

  static constexpr int V_ID_BITS = 21;
  static constexpr int W_BITS = 11;
  static constexpr int SHIFT_W_NEXT = V_ID_BITS;
  static constexpr int SHIFT_W = SHIFT_W_NEXT + W_BITS;
  static constexpr int SHIFT_LEVEL = SHIFT_W + W_BITS;

  uint64_t pack_key(int level, int w, int w_next, int v_id) const {
    return (static_cast<uint64_t>(level) << SHIFT_LEVEL) |
           (static_cast<uint64_t>(w) << SHIFT_W) |
           (static_cast<uint64_t>(w_next) << SHIFT_W_NEXT) |
           static_cast<uint64_t>(v_id);
  }

  void precompute_popcount() {
    popcnt.resize(B);
    for (int i = 0; i < B; ++i) {
      popcnt[i] = __builtin_popcount(static_cast<unsigned int>(i));
    }
  }

  void precompute_s_base() {
    s_base.resize(B + 1);
    s_base[0] = 0;
    for (int i = 0; i < B; ++i) {
      s_base[i + 1] = s_base[i] + popcnt[i];
    }
  }

  void precompute_tmax_in() {
    tmax_in.assign(n + 1, vector<uint16_t>(B, 0));
    for (int w = 0; w <= n; ++w) {
      for (int i = 0; i < B; ++i) {
        int maxStep = min(n, B - 1 - i);
        int lo = 0;
        int hi = maxStep;
        while (lo < hi) {
          int mid = (lo + hi + 1) / 2;
          int cost = (s_base[i + mid + 1] - s_base[i + 1]) + w * mid;
          if (cost <= n) {
            lo = mid;
          } else {
            hi = mid - 1;
          }
        }
        tmax_in[w][i] = static_cast<uint16_t>(lo);
      }
    }
  }

  void precompute_max_t_out() {
    max_t_out.assign(n + 2, vector<int16_t>(n + 1, -1));
    for (int w_next = 0; w_next <= n + 1; ++w_next) {
      int t = -1;
      for (int budget = 0; budget <= n; ++budget) {
        while (t + 1 < B) {
          int cand = t + 1;
          int cost = s_base[cand + 1] + w_next * (cand + 1);
          if (cost <= budget) {
            t = cand;
          } else {
            break;
          }
        }
        max_t_out[w_next][budget] = static_cast<int16_t>(t);
      }
    }
  }

  bool get_bit(const vector<uint64_t> &bits, int idx) const {
    return (bits[idx >> 6] >> (idx & 63)) & 1ULL;
  }

  void set_bit(vector<uint64_t> &bits, int idx) {
    bits[idx >> 6] |= (1ULL << (idx & 63));
  }

  int intern(const vector<uint64_t> &bits) {
    auto it = pool_map.find(bits);
    if (it != pool_map.end()) {
      return it->second;
    }
    int id = static_cast<int>(pool.size());
    pool.push_back(bits);
    pool_map.emplace(pool.back(), id);
    return id;
  }

  void bit_reset(int L) { fill(bit.begin(), bit.begin() + (L + 1), 0); }

  void bit_add(int idx, int delta, int L) {
    for (int i = idx + 1; i <= L; i += i & -i) {
      bit[i] += delta;
    }
  }

  int bit_sum(int idx) const {
    int res = 0;
    for (int i = idx + 1; i > 0; i -= i & -i) {
      res += bit[i];
    }
    return res;
  }

  int bit_range_sum(int l, int r) const {
    if (r < l) {
      return 0;
    }
    return bit_sum(r) - (l > 0 ? bit_sum(l - 1) : 0);
  }

  int compute_block(int w, int w_next, const vector<uint64_t> &v_out,
                    bool partial) {
    int L = partial ? (B - 1) : B;

    if (!partial) {
      uint8_t prev = 0;
      for (int i = 0; i < B; ++i) {
        bool lose = !get_bit(v_out, i);
        prev = static_cast<uint8_t>(prev | (lose ? 1 : 0));
        prefixLose[i] = prev;
      }
    }

    bit_reset(L);
    fill(status.begin(), status.end(), 0);

    for (int i = L - 1; i >= 0; --i) {
      bool win = false;

      int t_in = tmax_in[w][i];
      int max_in = i + t_in;
      if (max_in > L - 1) {
        max_in = L - 1;
      }
      if (max_in > i) {
        if (bit_range_sum(i + 1, max_in) > 0) {
          win = true;
        }
      }

      if (!win && !partial) {
        int cost_to_boundary = base_cost_to_end[i] + w * (L - 1 - i);
        if (cost_to_boundary <= n) {
          int budget = n - cost_to_boundary;
          int t_out = max_t_out[w_next][budget];
          if (t_out >= 0 && prefixLose[t_out]) {
            win = true;
          }
        }
      }

      status[i] = static_cast<uint8_t>(win ? 1 : 0);
      if (!win) {
        bit_add(i, 1, L);
      }
    }

    vector<uint64_t> bits(wordCount, 0ULL);
    for (int i = 0; i < L; ++i) {
      if (status[i]) {
        set_bit(bits, i);
      }
    }
    if (partial) {
      for (int i = L; i < B; ++i) {
        set_bit(bits, i);
      }
    }
    return intern(bits);
  }

  int solve_full(int lvl, int w, int w_next, int v_out_id) {
    if (w > n) {
      return 0;
    }
    uint64_t key = pack_key(lvl, w, w_next, v_out_id);
    auto it = memo_full.find(key);
    if (it != memo_full.end()) {
      return it->second;
    }

    int result;
    if (lvl == 0) {
      result = compute_block(w, w_next, pool[v_out_id], false);
    } else {
      int mid = solve_full(lvl - 1, w + 1, w_next, v_out_id);
      result = solve_full(lvl - 1, w, w + 1, mid);
    }
    memo_full.emplace(key, result);
    return result;
  }

  int solve_partial(int lvl, int w, int w_next, int v_out_id) {
    if (w > n) {
      return 0;
    }
    uint64_t key = pack_key(lvl, w, w_next, v_out_id);
    auto it = memo_partial.find(key);
    if (it != memo_partial.end()) {
      return it->second;
    }

    int result;
    if (lvl == 0) {
      result = compute_block(w, w_next, pool[v_out_id], true);
    } else {
      int mid = solve_partial(lvl - 1, w + 1, w_next, v_out_id);
      result = solve_full(lvl - 1, w, w + 1, mid);
    }
    memo_partial.emplace(key, result);
    return result;
  }
};

static bool validate_small() {
  struct Check {
    int n;
    int expected;
  };
  vector<Check> checks = {{2, 2}, {7, 1}, {20, 4}};
  for (const auto &c : checks) {
    Solver solver(c.n);
    long long m = solver.compute_M();
    if (m != c.expected) {
      cerr << "Validation failed: M(" << c.n << ") = " << m << ", expected "
           << c.expected << "\n";
      return false;
    }
  }

  long long sum20 = 0;
  for (int n = 1; n <= 20; ++n) {
    Solver solver(n);
    long long m = solver.compute_M();
    sum20 += m * m * m;
  }
  if (sum20 != 8150) {
    cerr << "Validation failed: sum_{n<=20} M(n)^3 = " << sum20
         << ", expected 8150\n";
    return false;
  }
  cout << "Validation passed (M(2)=2, M(7)=1, M(20)=4, sum20=8150).\n";
  return true;
}

int main() {
  if (!validate_small()) {
    return 1;
  }

  const int MAX_N = 1000;
  atomic<long long> total_sum{0};
  atomic<int> next_n{1};
  atomic<int> processed{0};
  mutex io_mutex;

  int num_threads = static_cast<int>(thread::hardware_concurrency());
  if (num_threads <= 0) {
    num_threads = 4;
  }

  num_threads = max(num_threads, 4);

  auto worker = [&]() {
    while (true) {
      int n = next_n.fetch_add(1);
      if (n > MAX_N) {
        break;
      }
      Solver solver(n);
      long long m = solver.compute_M();
      total_sum.fetch_add(m * m * m, memory_order_relaxed);

      int done = processed.fetch_add(1) + 1;
      if (done % 50 == 0) {
        lock_guard<mutex> lock(io_mutex);
        cout << "Processed " << done << "/" << MAX_N << "\n";
      }
    }
  };

  vector<thread> threads;
  threads.reserve(num_threads);
  for (int i = 0; i < num_threads; ++i) {
    threads.emplace_back(worker);
  }
  for (auto &t : threads) {
    t.join();
  }

  cout << "Final Result: " << total_sum.load() << "\n";
  return 0;
}
