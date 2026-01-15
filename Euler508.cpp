#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <future>
#include <iostream>
#include <thread>
#include <unordered_map>

using namespace std;

namespace {

constexpr int64_t MOD = 1000000007LL;
constexpr int64_t BRUTE_LIMIT = 1024;

struct Key {
    int64_t a = 0;
    int64_t b = 0;
    int64_t c = 0;
    int64_t d = 0;

    bool operator==(const Key& other) const {
        return a == other.a && b == other.b && c == other.c && d == other.d;
    }
};

struct KeyHash {
    size_t operator()(const Key& k) const {
        size_t h = std::hash<int64_t>{}(k.a);
        h ^= std::hash<int64_t>{}(k.b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<int64_t>{}(k.c) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<int64_t>{}(k.d) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

int64_t floor_div(int64_t a, int64_t b) {
    int64_t q = a / b;
    int64_t r = a % b;
    if (r != 0 && a < 0) --q;
    return q;
}

int64_t ceil_div(int64_t a, int64_t b) {
    return -floor_div(-a, b);
}

void count_even_odd(int64_t l, int64_t r, int64_t& even, int64_t& odd) {
    if (l > r) {
        even = 0;
        odd = 0;
        return;
    }
    int64_t first_even = (l % 2 == 0) ? l : l + 1;
    int64_t last_even = (r % 2 == 0) ? r : r - 1;
    if (first_even > last_even) {
        even = 0;
    } else {
        even = (last_even - first_even) / 2 + 1;
    }
    int64_t total = r - l + 1;
    odd = total - even;
}

int64_t mod_add(int64_t a, int64_t b) {
    a += b;
    if (a >= MOD) a -= MOD;
    return a;
}

int64_t ones_in_base(int64_t a, int64_t b) {
    int64_t count = 0;
    while (a != 0 || b != 0) {
        int64_t r = (a - b) & 1LL;
        count += r;
        int64_t na = (b - a + r) / 2;
        int64_t nb = (r - a - b) / 2;
        a = na;
        b = nb;
    }
    return count;
}

class Solver {
  public:
    int64_t compute_square(int64_t L) { return sum_b(-L, L, -L, L); }
    int64_t compute_a(int64_t umin, int64_t umax, int64_t vmin, int64_t vmax) {
        return sum_a(umin, umax, vmin, vmax);
    }

  private:
    unordered_map<Key, int64_t, KeyHash> memo_b_;
    unordered_map<Key, int64_t, KeyHash> memo_a_;

    int64_t brute_b(int64_t amin, int64_t amax, int64_t bmin, int64_t bmax) {
        int64_t total = 0;
        for (int64_t a = amin; a <= amax; ++a) {
            for (int64_t b = bmin; b <= bmax; ++b) {
                total += ones_in_base(a, b);
            }
        }
        return total % MOD;
    }

    int64_t sum_b(int64_t amin, int64_t amax, int64_t bmin, int64_t bmax) {
        if (amin > amax || bmin > bmax) return 0;
        Key key{amin, amax, bmin, bmax};
        auto it = memo_b_.find(key);
        if (it != memo_b_.end()) return it->second;

        __int128 area = static_cast<__int128>(amax - amin + 1) * static_cast<__int128>(bmax - bmin + 1);
        if (area <= BRUTE_LIMIT) {
            int64_t res = brute_b(amin, amax, bmin, bmax);
            memo_b_[key] = res;
            return res;
        }

        int64_t a_even, a_odd, b_even, b_odd;
        count_even_odd(amin, amax, a_even, a_odd);
        count_even_odd(bmin, bmax, b_even, b_odd);

        __int128 odd_pairs = static_cast<__int128>(a_even) * b_odd + static_cast<__int128>(a_odd) * b_even;
        int64_t res = static_cast<int64_t>(odd_pairs % MOD);

        int64_t u0_min = -amax;
        int64_t u0_max = -amin;
        int64_t u1_min = 1 - amax;
        int64_t u1_max = 1 - amin;

        res = mod_add(res, sum_a(u0_min, u0_max, bmin, bmax));
        res = mod_add(res, sum_a(u1_min, u1_max, bmin, bmax));

        memo_b_[key] = res;
        return res;
    }

    int64_t sum_a(int64_t umin, int64_t umax, int64_t vmin, int64_t vmax) {
        if (umin > umax || vmin > vmax) return 0;
        Key key{umin, umax, vmin, vmax};
        auto it = memo_a_.find(key);
        if (it != memo_a_.end()) return it->second;

        int64_t u_even, u_odd, v_even, v_odd;
        count_even_odd(umin, umax, u_even, u_odd);
        count_even_odd(vmin, vmax, v_even, v_odd);

        __int128 odd_pairs = static_cast<__int128>(u_odd) * v_odd;
        int64_t res = static_cast<int64_t>(odd_pairs % MOD);

        if (u_even && v_even) {
            int64_t x_min = ceil_div(umin, 2);
            int64_t x_max = floor_div(umax, 2);
            int64_t y_min = ceil_div(vmin, 2);
            int64_t y_max = floor_div(vmax, 2);
            if (x_min <= x_max && y_min <= y_max) {
                int64_t amin = -y_max;
                int64_t amax = -y_min;
                int64_t bmin = -x_max;
                int64_t bmax = -x_min;
                res = mod_add(res, sum_b(amin, amax, bmin, bmax));
            }
        }

        if (u_odd && v_odd) {
            int64_t x_min = ceil_div(umin - 1, 2);
            int64_t x_max = floor_div(umax - 1, 2);
            int64_t y_min = ceil_div(vmin - 1, 2);
            int64_t y_max = floor_div(vmax - 1, 2);
            if (x_min <= x_max && y_min <= y_max) {
                int64_t amin = -y_max;
                int64_t amax = -y_min;
                int64_t bmin = -x_max;
                int64_t bmax = -x_min;
                res = mod_add(res, sum_b(amin, amax, bmin, bmax));
            }
        }

        memo_a_[key] = res;
        return res;
    }
};

int64_t brute_square(int64_t L) {
    int64_t total = 0;
    for (int64_t a = -L; a <= L; ++a) {
        for (int64_t b = -L; b <= L; ++b) {
            total += ones_in_base(a, b);
        }
    }
    return total % MOD;
}

int64_t compute_parallel(int64_t L, unsigned threads) {
    if (threads <= 1) {
        Solver solver;
        return solver.compute_square(L);
    }

    int64_t amin = -L;
    int64_t amax = L;
    int64_t bmin = -L;
    int64_t bmax = L;

    int64_t a_even, a_odd, b_even, b_odd;
    count_even_odd(amin, amax, a_even, a_odd);
    count_even_odd(bmin, bmax, b_even, b_odd);

    __int128 odd_pairs = static_cast<__int128>(a_even) * b_odd + static_cast<__int128>(a_odd) * b_even;
    int64_t base = static_cast<int64_t>(odd_pairs % MOD);

    int64_t u0_min = -amax;
    int64_t u0_max = -amin;
    int64_t u1_min = 1 - amax;
    int64_t u1_max = 1 - amin;

    Solver solver0;
    Solver solver1;

    auto fut0 = async(launch::async, [&]() { return solver0.compute_a(u0_min, u0_max, bmin, bmax); });
    auto fut1 = async(launch::async, [&]() { return solver1.compute_a(u1_min, u1_max, bmin, bmax); });

    int64_t res = base;
    res = mod_add(res, fut0.get());
    res = mod_add(res, fut1.get());
    return res;
}

bool validate() {
    if (ones_in_base(11, 24) != 9) {
        cerr << "Validation failed: f(11+24i) mismatch." << '\n';
        return false;
    }
    if (ones_in_base(24, -11) != 7) {
        cerr << "Validation failed: f(24-11i) mismatch." << '\n';
        return false;
    }
    if (ones_in_base(8, 0) != 3) {
        cerr << "Validation failed: f(8+0i) mismatch." << '\n';
        return false;
    }
    if (ones_in_base(-5, 0) != 5) {
        cerr << "Validation failed: f(-5+0i) mismatch." << '\n';
        return false;
    }

    if (brute_square(1) != 20) {
        cerr << "Validation failed: B(1) mismatch." << '\n';
        return false;
    }
    if (brute_square(2) != 75) {
        cerr << "Validation failed: B(2) mismatch." << '\n';
        return false;
    }

    Solver solver;
    int64_t brute20 = brute_square(20);
    int64_t fast20 = solver.compute_square(20);
    if (fast20 != brute20) {
        cerr << "Validation failed: B(20) recursion mismatch." << '\n';
        return false;
    }
    int64_t fast500 = solver.compute_square(500);
    if (fast500 != 10795060) {
        cerr << "Validation failed: B(500) mismatch." << '\n';
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (!validate()) return 1;

    int64_t L = 1000000000000000LL;
    if (argc > 1) L = max<int64_t>(0, atoll(argv[1]));
    unsigned threads = thread::hardware_concurrency();
    if (argc > 2) threads = max(1, atoi(argv[2]));

    unsigned use_threads = min<unsigned>(2, threads);
    int64_t result = compute_parallel(L, use_threads);
    cout << result % MOD << '\n';
    return 0;
}
