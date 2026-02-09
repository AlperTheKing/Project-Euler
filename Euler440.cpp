#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr int MOD = 987898789;

struct Options {
    int limit = 2000;
    int threads = 0;
    bool run_checkpoints = true;
};

struct Mat2 {
    int a00 = 1;
    int a01 = 0;
    int a10 = 0;
    int a11 = 1;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1 && options.threads >= 0;
}

int choose_thread_count(int requested, std::size_t work_items) {
    if (work_items <= 1U) {
        return 1;
    }
    int threads = requested;
    if (threads <= 0) {
        threads = static_cast<int>(std::thread::hardware_concurrency());
    }
    if (threads <= 0) {
        threads = 4;
    }
    if (threads > static_cast<int>(work_items)) {
        threads = static_cast<int>(work_items);
    }
    if (threads < 1) {
        threads = 1;
    }
    return threads;
}

int mod_pow(int base, u64 exp, int mod) {
    u64 result = 1ULL;
    u64 cur = static_cast<u64>((base % mod + mod) % mod);
    u64 e = exp;
    while (e > 0ULL) {
        if (e & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * cur) % static_cast<u64>(mod));
        }
        cur = static_cast<u64>((static_cast<u128>(cur) * cur) % static_cast<u64>(mod));
        e >>= 1ULL;
    }
    return static_cast<int>(result);
}

Mat2 mul(const Mat2& x, const Mat2& y) {
    Mat2 r;
    r.a00 = static_cast<int>((static_cast<u128>(x.a00) * y.a00 + static_cast<u128>(x.a01) * y.a10) % MOD);
    r.a01 = static_cast<int>((static_cast<u128>(x.a00) * y.a01 + static_cast<u128>(x.a01) * y.a11) % MOD);
    r.a10 = static_cast<int>((static_cast<u128>(x.a10) * y.a00 + static_cast<u128>(x.a11) * y.a10) % MOD);
    r.a11 = static_cast<int>((static_cast<u128>(x.a10) * y.a01 + static_cast<u128>(x.a11) * y.a11) % MOD);
    return r;
}

Mat2 mat_pow(Mat2 base, u64 exp) {
    Mat2 result;
    u64 e = exp;
    while (e > 0ULL) {
        if (e & 1ULL) {
            result = mul(result, base);
        }
        e >>= 1ULL;
        if (e > 0ULL) {
            base = mul(base, base);
        }
    }
    return result;
}

std::vector<std::pair<u64, int>> factorize(u64 n) {
    std::vector<std::pair<u64, int>> fac;
    u64 x = n;
    for (u64 p = 2; p * p <= x; ++p) {
        if (x % p != 0ULL) {
            continue;
        }
        int e = 0;
        while (x % p == 0ULL) {
            x /= p;
            ++e;
        }
        fac.push_back({p, e});
    }
    if (x > 1ULL) {
        fac.push_back({x, 1});
    }
    return fac;
}

u64 sequence_period() {
    // T(n) recurrence: T_n = 10*T_{n-1} + T_{n-2} (mod MOD).
    // Companion matrix A = [[10,1],[1,0]].
    const Mat2 A{10, 1, 1, 0};
    const int leg = mod_pow(104 % MOD, static_cast<u64>((MOD - 1) / 2), MOD);
    u64 candidate = (leg == 1) ? static_cast<u64>(MOD - 1) : static_cast<u64>(MOD + 1);

    u64 period = candidate;
    const auto fac = factorize(candidate);
    for (const auto& [q, _e] : fac) {
        while (period % q == 0ULL) {
            const u64 trial = period / q;
            const Mat2 p = mat_pow(A, trial);
            if (p.a00 == 1 && p.a01 == 0 && p.a10 == 0 && p.a11 == 1) {
                period = trial;
            } else {
                break;
            }
        }
    }
    return period;
}

std::vector<Mat2> build_powers_for_T() {
    std::vector<Mat2> pw;
    pw.reserve(64);
    pw.push_back(Mat2{10, 1, 1, 0});
    for (int i = 1; i < 64; ++i) {
        pw.push_back(mul(pw.back(), pw.back()));
    }
    return pw;
}

int T_mod_from_index(u64 n, const std::vector<Mat2>& pw) {
    if (n == 0ULL) {
        return 1;  // T(0)
    }

    int v0 = 10;  // T(1)
    int v1 = 1;   // T(0)
    u64 e = n - 1ULL;
    int bit = 0;
    while (e > 0ULL) {
        if (e & 1ULL) {
            const Mat2& m = pw[static_cast<std::size_t>(bit)];
            const int nv0 =
                static_cast<int>((static_cast<u128>(m.a00) * v0 + static_cast<u128>(m.a01) * v1) % MOD);
            const int nv1 =
                static_cast<int>((static_cast<u128>(m.a10) * v0 + static_cast<u128>(m.a11) * v1) % MOD);
            v0 = nv0;
            v1 = nv1;
        }
        e >>= 1ULL;
        ++bit;
    }
    return v0;
}

std::vector<i64> build_counts_same(int l) {
    std::vector<i64> cnt(static_cast<std::size_t>(l + 1), 0);
    std::vector<int> v2(static_cast<std::size_t>(l + 1), 0);
    for (int x = 1; x <= l; ++x) {
        v2[static_cast<std::size_t>(x)] = __builtin_ctz(static_cast<unsigned>(x));
    }

    for (int a = 1; a <= l; ++a) {
        const int ta = v2[static_cast<std::size_t>(a)];
        for (int b = 1; b <= l; ++b) {
            if (v2[static_cast<std::size_t>(b)] != ta) {
                continue;
            }
            const int g = std::gcd(a, b);
            ++cnt[static_cast<std::size_t>(g)];
        }
    }
    return cnt;
}

int solve_mod(int l, int requested_threads) {
    const u64 period = sequence_period();
    const std::vector<Mat2> pw = build_powers_for_T();
    const std::vector<i64> counts_same = build_counts_same(l);

    i64 same_total = 0;
    for (int d = 1; d <= l; ++d) {
        same_total += counts_same[static_cast<std::size_t>(d)];
    }
    const i64 total_pairs = static_cast<i64>(l) * static_cast<i64>(l);
    const i64 count_diff = total_pairs - same_total;

    std::vector<int> counts_same_mod(static_cast<std::size_t>(l + 1), 0);
    for (int d = 1; d <= l; ++d) {
        counts_same_mod[static_cast<std::size_t>(d)] =
            static_cast<int>(counts_same[static_cast<std::size_t>(d)] % MOD);
    }

    const int thread_count = choose_thread_count(requested_threads, static_cast<std::size_t>(l));
    std::vector<int> partial(static_cast<std::size_t>(thread_count), 0);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));

    for (int t = 0; t < thread_count; ++t) {
        const int begin = l * t / thread_count + 1;
        const int end = l * (t + 1) / thread_count;

        workers.emplace_back([&, begin, end, t]() {
            int local = 0;

            for (int c = begin; c <= end; ++c) {
                const int base = (c & 1) ? 10 : 1;
                int sum_c = static_cast<int>((static_cast<u128>(count_diff % MOD) * base) % MOD);

                u64 p = 1ULL;
                for (int d = 1; d <= l; ++d) {
                    p = static_cast<u64>((static_cast<u128>(p) * static_cast<u64>(c)) % period);
                    const int tval = T_mod_from_index(p, pw);
                    sum_c = static_cast<int>((sum_c +
                                              static_cast<u128>(counts_same_mod[static_cast<std::size_t>(d)]) * tval) %
                                             MOD);
                }

                local += sum_c;
                if (local >= MOD) {
                    local -= MOD;
                }
            }

            partial[static_cast<std::size_t>(t)] = local;
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }

    int ans = 0;
    for (int x : partial) {
        ans += x;
        if (ans >= MOD) {
            ans -= MOD;
        }
    }
    return ans;
}

bool run_checkpoints(int requested_threads) {
    if (solve_mod(2, requested_threads) != 10444) {
        std::cerr << "Checkpoint failed: S(2)\n";
        return false;
    }
    if (solve_mod(3, requested_threads) != 862246950) {
        std::cerr << "Checkpoint failed: S(3) mod 987898789\n";
        return false;
    }
    if (solve_mod(4, requested_threads) != 670616280) {
        std::cerr << "Checkpoint failed: S(4) mod 987898789\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints(options.threads)) {
        return 2;
    }

    std::cout << solve_mod(options.limit, options.threads) << '\n';
    return 0;
}
