#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

struct Options {
    int limit = 1000000;
    unsigned threads = std::thread::hardware_concurrency();
    bool run_checkpoints = true;
};

struct PrimeExp {
    int p;
    int e;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg, const std::string& prefix,
                                 unsigned& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    unsigned parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10U + static_cast<unsigned>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.limit < 3) {
        return false;
    }
    if (options.threads == 0) {
        options.threads = 1;
    }
    return true;
}

int iroot_floor(i64 x) {
    if (x <= 0) {
        return 0;
    }
    i64 r = static_cast<i64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r > x) {
        --r;
    }
    return static_cast<int>(r);
}

int iroot_ceil(i64 x) {
    if (x <= 0) {
        return 0;
    }
    int r = iroot_floor(x);
    if (static_cast<i64>(r) * r < x) {
        ++r;
    }
    return r;
}

int cube_root_floor(i64 x) {
    int lo = 0;
    int hi = static_cast<int>(std::cbrt(static_cast<long double>(x))) + 5;
    while (static_cast<i64>(hi) * hi * hi <= x) {
        ++hi;
    }
    while (lo + 1 < hi) {
        const int mid = lo + (hi - lo) / 2;
        const i64 cube = static_cast<i64>(mid) * mid * mid;
        if (cube <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

std::vector<int> build_spf(int n) {
    std::vector<int> spf(n + 1);
    for (int i = 0; i <= n; ++i) {
        spf[i] = i;
    }
    if (n >= 0) {
        spf[0] = 0;
    }
    if (n >= 1) {
        spf[1] = 1;
    }

    for (int i = 2; static_cast<i64>(i) * i <= n; ++i) {
        if (spf[i] != i) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            if (spf[j] == j) {
                spf[j] = i;
            }
        }
    }
    return spf;
}

std::vector<int> build_squarefree_parts(int n, const std::vector<int>& spf) {
    std::vector<int> sf(n + 1, 1);
    sf[0] = 0;
    for (int x = 1; x <= n; ++x) {
        int t = x;
        int value = 1;
        while (t > 1) {
            const int p = spf[t];
            int c = 0;
            while (t % p == 0) {
                t /= p;
                c ^= 1;
            }
            if (c != 0) {
                value *= p;
            }
        }
        sf[x] = value;
    }
    return sf;
}

std::vector<std::vector<PrimeExp>> build_factor_table(int n, const std::vector<int>& spf) {
    std::vector<std::vector<PrimeExp>> table(n + 1);
    for (int x = 2; x <= n; ++x) {
        int t = x;
        while (t > 1) {
            const int p = spf[t];
            int e = 0;
            while (t % p == 0) {
                t /= p;
                ++e;
            }
            table[x].push_back({p, e});
        }
    }
    return table;
}

bool multiply_power_limited(i64& value, int p, int e, int limit) {
    for (int i = 0; i < e; ++i) {
        value *= p;
        if (value > limit) {
            return false;
        }
    }
    return true;
}

u64 count_for_k(int limit, int k, const std::vector<int>& squarefree_u,
                const std::vector<std::vector<PrimeExp>>& factors) {
    const int L = limit / k;
    if (L < 2) {
        return 0;
    }

    u64 subtotal = 0;
    const int v_min = (k + 1) / 2;

    for (int v = v_min; v < k; ++v) {
        if (std::gcd(v, k) != 1) {
            continue;
        }

        const int u = k - v;
        const int sf_u = squarefree_u[u];
        const i64 W = static_cast<i64>(v) * sf_u;
        if (W > static_cast<i64>(L) * L) {
            continue;
        }

        std::vector<PrimeExp> wf = factors[v];
        if (sf_u != 1) {
            for (const PrimeExp pe : factors[sf_u]) {
                wf.push_back({pe.p, 1});
            }
        }
        std::sort(wf.begin(), wf.end(), [](const PrimeExp& a, const PrimeExp& b) {
            return a.p < b.p;
        });

        const int r_max = static_cast<int>((static_cast<i64>(u) * L) / (u + 2 * v));
        for (int r = 1; r <= r_max; ++r) {
            const std::vector<PrimeExp>& rf = factors[r];

            i64 s0 = 1;
            std::size_t i = 0;
            std::size_t j = 0;
            bool valid = true;

            while (i < wf.size() || j < rf.size()) {
                int p = 0;
                int e = 0;
                int a = 0;

                if (j == rf.size() || (i < wf.size() && wf[i].p < rf[j].p)) {
                    p = wf[i].p;
                    e = wf[i].e;
                    a = 0;
                    ++i;
                } else if (i == wf.size() || rf[j].p < wf[i].p) {
                    p = rf[j].p;
                    e = 0;
                    a = rf[j].e;
                    ++j;
                } else {
                    p = wf[i].p;
                    e = wf[i].e;
                    a = rf[j].e;
                    ++i;
                    ++j;
                }

                int b0 = 0;
                if (a < e) {
                    b0 = e - a;
                } else {
                    b0 = (a - e) & 1;
                }

                if (b0 > 0) {
                    if (!multiply_power_limited(s0, p, b0, L)) {
                        valid = false;
                        break;
                    }
                }
            }

            if (!valid) {
                continue;
            }

            const int s_min = static_cast<int>(
                (static_cast<i64>(u + 2 * v) * r + (u - 1)) / u);
            if (s_min > L) {
                continue;
            }

            const i64 ratio_min = (static_cast<i64>(s_min) + s0 - 1) / s0;
            int n_min = iroot_ceil(ratio_min);
            const int n_max = iroot_floor(L / s0);
            if (n_min > n_max) {
                continue;
            }

            if ((s0 & 1LL) == 0) {
                if ((r & 1) == 0) {
                    subtotal += static_cast<u64>(n_max - n_min + 1);
                }
                continue;
            }

            const int want_parity = r & 1;
            if ((n_min & 1) != want_parity) {
                ++n_min;
            }
            if (n_min <= n_max) {
                subtotal += static_cast<u64>((n_max - n_min) / 2 + 1);
            }
        }
    }

    return subtotal;
}

u64 solve_fast(int limit, unsigned threads) {
    const i64 n2_twice = 2LL * limit * limit;
    const int k_max = cube_root_floor(n2_twice) + 2;
    const int max_needed = std::max(k_max, limit / 6 + 5);

    const std::vector<int> spf = build_spf(max_needed);
    const std::vector<int> squarefree_u = build_squarefree_parts(k_max, spf);
    const std::vector<std::vector<PrimeExp>> factors = build_factor_table(max_needed, spf);

    if (threads <= 1 || k_max < 128) {
        u64 total = 0;
        for (int k = 2; k <= k_max; ++k) {
            total += count_for_k(limit, k, squarefree_u, factors);
        }
        return total;
    }

    const unsigned use_threads = std::min<unsigned>(threads, static_cast<unsigned>(k_max - 1));
    std::atomic<int> next_k(2);
    std::vector<std::thread> pool;
    std::vector<u64> partial(use_threads, 0);
    pool.reserve(use_threads);

    for (unsigned t = 0; t < use_threads; ++t) {
        pool.emplace_back([&, t]() {
            u64 local = 0;
            while (true) {
                const int k = next_k.fetch_add(1, std::memory_order_relaxed);
                if (k > k_max) {
                    break;
                }
                local += count_for_k(limit, k, squarefree_u, factors);
            }
            partial[t] = local;
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    u64 total = 0;
    for (u64 part : partial) {
        total += part;
    }
    return total;
}

u64 solve_bruteforce(int limit) {
    u64 total = 0;
    for (int a = 1; a <= limit / 3; ++a) {
        for (int b = a; b <= (limit - a) / 2; ++b) {
            const int s = a + b;
            const int c_max = std::min(s - 1, limit - s);
            for (int c = b; c <= c_max; ++c) {
                const i64 num = static_cast<i64>(a) * a * a * (s - c) * (s + c);
                const i64 den = static_cast<i64>(b) * s * s;
                if (num % den != 0) {
                    continue;
                }
                const i64 q = num / den;
                const i64 rt = static_cast<i64>(std::sqrt(static_cast<long double>(q)));
                if (rt * rt == q || (rt + 1) * (rt + 1) == q) {
                    ++total;
                }
            }
        }
    }
    return total;
}

bool run_checkpoints() {
    for (int limit : {50, 80, 120}) {
        const u64 fast = solve_fast(limit, 1);
        const u64 brute = solve_bruteforce(limit);
        if (fast != brute) {
            std::cerr << "Bruteforce checkpoint failed for limit=" << limit << ": fast=" << fast
                      << ", brute=" << brute << '\n';
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const u64 answer = solve_fast(options.limit, options.threads);
    std::cout << answer << '\n';
    return 0;
}
