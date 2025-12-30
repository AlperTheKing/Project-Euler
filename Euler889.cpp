#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

namespace {

constexpr int64_t MOD = 1'000'062'031LL;

int64_t mod_mul(int64_t a, int64_t b) {
    return static_cast<int64_t>((__int128)a * b % MOD);
}

int64_t mod_pow(int64_t base, uint64_t exp) {
    int64_t result = 1 % MOD;
    int64_t cur = base % MOD;
    while (exp > 0) {
        if (exp & 1) result = mod_mul(result, cur);
        cur = mod_mul(cur, cur);
        exp >>= 1;
    }
    return result;
}

int64_t mod_add(int64_t a, int64_t b) {
    int64_t v = a + b;
    if (v >= MOD) v -= MOD;
    return v;
}

int64_t mod_sub(int64_t a, int64_t b) {
    int64_t v = a - b;
    if (v < 0) v += MOD;
    return v;
}

int64_t brute_mod(int64_t k, int64_t t, int r) {
    cpp_int A = (cpp_int(1) << t) + 1;
    cpp_int Apow = 1;
    for (int i = 0; i < r; ++i) {
        Apow *= A;
    }
    cpp_int B = (cpp_int(1) << k) + 1;
    cpp_int total = 0;
    for (int64_t n = 0; n < k; ++n) {
        cpp_int rn = (Apow << n) % B;
        cpp_int dn = rn;
        cpp_int comp = B - rn;
        if (comp < dn) dn = comp;
        total += dn << (k - n);
    }
    int64_t result = (total % MOD).convert_to<long long>();
    if (result < 0) result += MOD;
    return result;
}

struct FastContext {
    int r = 0;
    int64_t k = 0;
    int64_t t = 0;
    vector<uint64_t> binom;
    vector<int64_t> binom_mod;
    vector<int64_t> pow2_ti;
    int64_t pow2_k = 0;
    int64_t pow2_k_inv = 0;
    int64_t b_mod = 0;
    vector<int64_t> prefix_k;
    vector<int64_t> suffix;
};

struct Task {
    int64_t n = 0;
    int s = 0;
};

int64_t exception_contrib(const FastContext& ctx, int64_t n, int s) {
    int64_t pow2_n = mod_pow(2, static_cast<uint64_t>(n));
    int64_t S_mod = 0;
    for (int i = 0; i <= ctx.r; ++i) {
        int64_t base = mod_mul(ctx.binom_mod[i], mod_mul(pow2_n, ctx.pow2_ti[i]));
        if (i < s) {
            S_mod = mod_add(S_mod, base);
        } else {
            int64_t term = mod_mul(base, ctx.pow2_k_inv);
            S_mod = mod_sub(S_mod, term);
        }
    }

    int idx = (s <= ctx.r) ? (s - 1) : ctx.r;
    int64_t e_max = n + ctx.t * static_cast<int64_t>(idx);
    int64_t shift = ctx.k - e_max;
    uint64_t C_top = ctx.binom[idx];

    uint64_t Q = 0;
    if (shift >= 0 && shift < 63) {
        Q = C_top >> shift;
    }
    int64_t r_mod = mod_sub(S_mod, mod_mul(static_cast<int64_t>(Q % MOD), ctx.b_mod));

    bool greater = false;
    if (shift < 63) {
        uint64_t low_mask = (static_cast<uint64_t>(1) << shift) - 1;
        uint64_t low = C_top & low_mask;
        uint64_t half = static_cast<uint64_t>(1) << (shift - 1);
        if (low < half) {
            greater = false;
        } else if (low > half) {
            greater = true;
        } else {
            greater = (s >= 2);
        }
    }

    int64_t d_mod = greater ? mod_sub(ctx.b_mod, r_mod) : r_mod;
    int64_t pow2_kn = mod_pow(2, static_cast<uint64_t>(ctx.k - n));
    return mod_mul(d_mod, pow2_kn);
}

int64_t compute_fast(int64_t k, int64_t t, int r) {
    FastContext ctx;
    ctx.r = r;
    ctx.k = k;
    ctx.t = t;
    ctx.binom.assign(r + 1, 0);
    ctx.binom_mod.assign(r + 1, 0);
    ctx.pow2_ti.assign(r + 1, 0);
    ctx.prefix_k.assign(r + 2, 0);
    ctx.suffix.assign(r + 2, 0);

    ctx.binom[0] = 1;
    for (int i = 1; i <= r; ++i) {
        __int128 num = static_cast<__int128>(ctx.binom[i - 1]) * (r - i + 1);
        ctx.binom[i] = static_cast<uint64_t>(num / i);
    }
    for (int i = 0; i <= r; ++i) {
        ctx.binom_mod[i] = static_cast<int64_t>(ctx.binom[i] % MOD);
    }

    int64_t pow2_t = mod_pow(2, static_cast<uint64_t>(t));
    ctx.pow2_ti[0] = 1;
    for (int i = 1; i <= r; ++i) {
        ctx.pow2_ti[i] = mod_mul(ctx.pow2_ti[i - 1], pow2_t);
    }

    ctx.pow2_k = mod_pow(2, static_cast<uint64_t>(k));
    ctx.pow2_k_inv = mod_pow(ctx.pow2_k, MOD - 2);
    ctx.b_mod = mod_add(ctx.pow2_k, 1);

    vector<int64_t> term(r + 1, 0);
    vector<int64_t> term_k(r + 1, 0);
    for (int i = 0; i <= r; ++i) {
        term[i] = mod_mul(ctx.binom_mod[i], ctx.pow2_ti[i]);
        term_k[i] = mod_mul(term[i], ctx.pow2_k);
    }

    for (int i = 0; i <= r; ++i) {
        ctx.prefix_k[i + 1] = mod_add(ctx.prefix_k[i], term_k[i]);
    }
    for (int i = r; i >= 0; --i) {
        ctx.suffix[i] = mod_add(ctx.suffix[i + 1], term[i]);
    }

    vector<Task> tasks;
    tasks.reserve(static_cast<size_t>((r + 1) * 62));
    int64_t total = 0;

    auto add_bulk = [&](int64_t bulk, int64_t constant) {
        if (bulk <= 0) return;
        int64_t bulk_mod = static_cast<int64_t>(bulk % MOD);
        total = mod_add(total, mod_mul(bulk_mod, constant));
    };

    for (int s = 1; s <= r; ++s) {
        __int128 start_raw = static_cast<__int128>(k) - static_cast<__int128>(s) * t;
        __int128 end_raw = static_cast<__int128>(k) - static_cast<__int128>(s - 1) * t - 1;
        if (end_raw < 0) continue;
        int64_t start = start_raw > 0 ? static_cast<int64_t>(start_raw) : 0;
        int64_t end = static_cast<int64_t>(end_raw);
        if (end > k - 1) end = k - 1;
        if (start > end) continue;

        int64_t length = end - start + 1;
        int64_t exc_start = max(start, end - 61);
        int64_t exc_count = (exc_start <= end) ? (end - exc_start + 1) : 0;
        int64_t bulk = length - exc_count;

        int64_t constant = mod_sub(ctx.prefix_k[s], ctx.suffix[s]);
        add_bulk(bulk, constant);

        for (int64_t n = exc_start; n <= end; ++n) {
            tasks.push_back({n, s});
        }
    }

    {
        int s = r + 1;
        __int128 end_raw = static_cast<__int128>(k) - static_cast<__int128>(r) * t - 1;
        if (end_raw >= 0) {
            int64_t start = 0;
            int64_t end = static_cast<int64_t>(end_raw);
            if (end > k - 1) end = k - 1;
            if (start <= end) {
                int64_t length = end - start + 1;
                int64_t exc_start = max(start, end - 61);
                int64_t exc_count = (exc_start <= end) ? (end - exc_start + 1) : 0;
                int64_t bulk = length - exc_count;

                int64_t constant = ctx.prefix_k[r + 1];
                add_bulk(bulk, constant);

                for (int64_t n = exc_start; n <= end; ++n) {
                    tasks.push_back({n, s});
                }
            }
        }
    }

    size_t task_count = tasks.size();
    if (task_count > 0) {
        int thread_count = static_cast<int>(thread::hardware_concurrency());
        if (thread_count <= 0) thread_count = 1;
        if (thread_count > static_cast<int>(task_count)) {
            thread_count = static_cast<int>(task_count);
        }
        vector<int64_t> partial(thread_count, 0);
        vector<thread> workers;
        workers.reserve(thread_count);

        for (int t_idx = 0; t_idx < thread_count; ++t_idx) {
            workers.emplace_back([&, t_idx]() {
                int64_t acc = 0;
                for (size_t i = t_idx; i < task_count; i += thread_count) {
                    acc += exception_contrib(ctx, tasks[i].n, tasks[i].s);
                    if (acc >= MOD) acc -= MOD;
                }
                partial[t_idx] = acc;
            });
        }
        for (auto& th : workers) th.join();
        for (int t_idx = 0; t_idx < thread_count; ++t_idx) {
            total = mod_add(total, partial[t_idx]);
        }
    }

    return total % MOD;
}

}  // namespace

int main() {
    if (brute_mod(3, 1, 1) != 42) {
        cerr << "Validation failure: F(3,1,1)\n";
        return 1;
    }
    if (brute_mod(13, 3, 3) != 23093880) {
        cerr << "Validation failure: F(13,3,3)\n";
        return 1;
    }
    if (brute_mod(103, 13, 6) != 878922518) {
        cerr << "Validation failure: F(103,13,6)\n";
        return 1;
    }

    const int64_t k = 1'000'000'000'000'000'000LL + 31;
    const int64_t t = 100'000'000'000'000LL + 31;
    const int r = 62;

    cout << compute_fast(k, t, r) << '\n';
    return 0;
}
