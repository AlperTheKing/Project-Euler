#include <algorithm>
#include <array>
#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using boost::multiprecision::cpp_int;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int MAXP = 20;

struct Triple {
    cpp_int a = 0;
    cpp_int b = 0;
    cpp_int c = 0;
    cpp_int sum = 0;
};

cpp_int cube_cpp(const u64 x) {
    cpp_int t = x;
    return t * t * t;
}

long double lower_bound_range(const long double l1, const long double l2, const long double l3, long double rem) {
    const long double t = l2 - l1;
    long double min_value;
    if (rem <= t) {
        min_value = l1 + rem;
    } else {
        rem -= t;
        const long double need = 2.0L * (l3 - l2);
        if (rem <= need) {
            min_value = l2 + rem / 2.0L;
        } else {
            min_value = l3;
        }
    }
    return l3 - min_value;
}

std::string to_string_cpp(const cpp_int& v) {
    return v.convert_to<std::string>();
}

struct Solver {
    std::vector<u64> primes;
    std::vector<int> exps;
    int pcount = 0;

    cpp_int n_value = 1;
    long double log_n = 0.0L;
    std::vector<long double> log_p;
    std::vector<long double> rem_log;

    std::vector<std::vector<u64>> pow_u64;
    std::vector<std::vector<u128>> pow_u128;
    std::vector<std::vector<cpp_int>> pow_cpp;

    std::vector<std::vector<std::array<int, 3>>> split_by_exp;
    std::vector<int> desc_idx;

    long long phase_a_nodes = 0;
    long long phase_a_cap = 3000000;
    long double best_log_range = 1e100L;
    Triple best;

    explicit Solver(std::vector<u64> primes_in, std::vector<int> exps_in)
        : primes(std::move(primes_in)), exps(std::move(exps_in)), pcount(static_cast<int>(primes.size())) {
        log_p.resize(pcount);
        rem_log.assign(pcount + 1, 0.0L);
        pow_u64.resize(pcount);
        pow_u128.resize(pcount);
        pow_cpp.resize(pcount);
        desc_idx.resize(pcount);
        std::iota(desc_idx.begin(), desc_idx.end(), 0);
        std::sort(desc_idx.begin(), desc_idx.end(), [&](const int i, const int j) { return primes[i] > primes[j]; });

        int max_e = 0;
        for (int i = 0; i < pcount; ++i) {
            log_p[i] = std::log(static_cast<long double>(primes[i]));
            log_n += static_cast<long double>(exps[i]) * log_p[i];
            n_value *= cpp_int(1);
            max_e = std::max(max_e, exps[i]);

            pow_u64[i].assign(static_cast<std::size_t>(exps[i] + 1), 1ULL);
            pow_u128[i].assign(static_cast<std::size_t>(exps[i] + 1), 1U);
            pow_cpp[i].assign(static_cast<std::size_t>(exps[i] + 1), cpp_int(1));
            for (int e = 1; e <= exps[i]; ++e) {
                pow_u64[i][static_cast<std::size_t>(e)] =
                    pow_u64[i][static_cast<std::size_t>(e - 1)] * primes[i];
                pow_u128[i][static_cast<std::size_t>(e)] =
                    pow_u128[i][static_cast<std::size_t>(e - 1)] * static_cast<u128>(primes[i]);
                pow_cpp[i][static_cast<std::size_t>(e)] =
                    pow_cpp[i][static_cast<std::size_t>(e - 1)] * cpp_int(primes[i]);
            }
            n_value *= pow_cpp[i][static_cast<std::size_t>(exps[i])];
        }

        for (int i = pcount - 1; i >= 0; --i) {
            rem_log[static_cast<std::size_t>(i)] =
                rem_log[static_cast<std::size_t>(i + 1)] + static_cast<long double>(exps[i]) * log_p[i];
        }

        split_by_exp.assign(static_cast<std::size_t>(max_e + 1), {});
        for (int e = 0; e <= max_e; ++e) {
            auto& v = split_by_exp[static_cast<std::size_t>(e)];
            for (int x = 0; x <= e; ++x) {
                for (int y = 0; y <= e - x; ++y) {
                    int z = e - x - y;
                    v.push_back({x, y, z});
                }
            }
            std::sort(v.begin(), v.end(), [&](const auto& a, const auto& b) {
                const int sa = std::max({a[0], a[1], a[2]}) - std::min({a[0], a[1], a[2]});
                const int sb = std::max({b[0], b[1], b[2]}) - std::min({b[0], b[1], b[2]});
                if (sa != sb) return sa < sb;
                const int ma = std::abs(3 * a[0] - e) + std::abs(3 * a[1] - e) + std::abs(3 * a[2] - e);
                const int mb = std::abs(3 * b[0] - e) + std::abs(3 * b[1] - e) + std::abs(3 * b[2] - e);
                return ma < mb;
            });
        }
    }

    static void sort_bins(std::array<long double, 3>& logs,
                          std::array<std::array<std::uint8_t, MAXP>, 3>& bins) {
        std::array<int, 3> ord = {0, 1, 2};
        std::sort(ord.begin(), ord.end(), [&](const int i, const int j) { return logs[i] < logs[j]; });
        std::array<long double, 3> nl;
        std::array<std::array<std::uint8_t, MAXP>, 3> nb;
        for (int i = 0; i < 3; ++i) {
            nl[i] = logs[ord[i]];
            nb[i] = bins[ord[i]];
        }
        logs = nl;
        bins = nb;
    }

    cpp_int value_from_bin(const std::array<std::uint8_t, MAXP>& bin_exp) const {
        cpp_int v = 1;
        for (int i = 0; i < pcount; ++i) {
            const int e = static_cast<int>(bin_exp[i]);
            if (e > 0) v *= pow_cpp[i][static_cast<std::size_t>(e)];
        }
        return v;
    }

    void update_best_from_bins(const std::array<long double, 3>& logs,
                               const std::array<std::array<std::uint8_t, MAXP>, 3>& bins) {
        const long double range = logs[2] - logs[0];
        if (range > best_log_range + 1e-18L) return;

        std::array<cpp_int, 3> v = {value_from_bin(bins[0]), value_from_bin(bins[1]), value_from_bin(bins[2])};
        std::sort(v.begin(), v.end());
        const cpp_int sum = v[0] + v[1] + v[2];

        if (range < best_log_range - 1e-18L || best.sum == 0 || sum < best.sum) {
            best_log_range = range;
            best.a = v[0];
            best.b = v[1];
            best.c = v[2];
            best.sum = sum;
        }
    }

    void greedy_seed() {
        std::array<long double, 3> logs = {0.0L, 0.0L, 0.0L};
        std::array<std::array<std::uint8_t, MAXP>, 3> bins{};
        for (auto& row : bins) row.fill(0U);

        std::vector<int> idx(pcount);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](const int i, const int j) { return primes[i] > primes[j]; });

        for (const int pi : idx) {
            for (int c = 0; c < exps[pi]; ++c) {
                int target = 0;
                if (logs[1] < logs[target]) target = 1;
                if (logs[2] < logs[target]) target = 2;
                logs[target] += log_p[pi];
                ++bins[target][pi];
            }
        }
        sort_bins(logs, bins);
        best_log_range = logs[2] - logs[0];
        update_best_from_bins(logs, bins);
    }

    void phase_a_dfs(int idx, const std::array<long double, 3>& logs,
                     const std::array<std::array<std::uint8_t, MAXP>, 3>& bins) {
        if (phase_a_nodes++ >= phase_a_cap) return;
        if (idx == pcount) {
            update_best_from_bins(logs, bins);
            return;
        }
        if (lower_bound_range(logs[0], logs[1], logs[2], rem_log[static_cast<std::size_t>(idx)]) >=
            best_log_range - 1e-18L) {
            return;
        }

        const int e = exps[idx];
        const auto& splits = split_by_exp[static_cast<std::size_t>(e)];
        for (const auto& s : splits) {
            std::array<long double, 3> nl = logs;
            std::array<std::array<std::uint8_t, MAXP>, 3> nb = bins;

            nl[0] += static_cast<long double>(s[0]) * log_p[idx];
            nl[1] += static_cast<long double>(s[1]) * log_p[idx];
            nl[2] += static_cast<long double>(s[2]) * log_p[idx];

            nb[0][idx] = static_cast<std::uint8_t>(nb[0][idx] + s[0]);
            nb[1][idx] = static_cast<std::uint8_t>(nb[1][idx] + s[1]);
            nb[2][idx] = static_cast<std::uint8_t>(nb[2][idx] + s[2]);

            sort_bins(nl, nb);
            phase_a_dfs(idx + 1, nl, nb);
        }
    }

    u64 cube_root_floor() const {
        long double approx = std::exp(log_n / 3.0L);
        if (approx < 1.0L) approx = 1.0L;
        u64 x = static_cast<u64>(approx);
        while (cube_cpp(x + 1) <= n_value) ++x;
        while (cube_cpp(x) > n_value) --x;
        return x;
    }

    u64 sqrt_floor_cpp(const cpp_int& x) const {
        long double approx = std::sqrt(x.convert_to<long double>());
        if (approx < 0.0L) approx = 0.0L;
        u64 r = static_cast<u64>(approx);
        while (cpp_int(r + 1U) * cpp_int(r + 1U) <= x) ++r;
        while (cpp_int(r) * cpp_int(r) > x) --r;
        return r;
    }

    u64 lower_a_bound(const u64 high) const {
        const cpp_int rhs = n_value * best.a * best.a;
        const cpp_int c2 = best.c * best.c;
        u64 lo = 1;
        u64 hi = high;
        while (lo < hi) {
            const u64 mid = lo + (hi - lo) / 2;
            const cpp_int lhs = c2 * cube_cpp(mid);
            if (lhs >= rhs) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        return lo;
    }

    u64 max_divisor_leq(const std::array<int, MAXP>& residual, const u64 limit) const {
        std::array<int, MAXP> e_desc{};
        for (int i = 0; i < pcount; ++i) {
            e_desc[i] = residual[desc_idx[i]];
        }

        std::array<u128, MAXP + 1> rem_max{};
        rem_max[pcount] = 1;
        for (int i = pcount - 1; i >= 0; --i) {
            const int pi = desc_idx[i];
            rem_max[i] = rem_max[i + 1] * pow_u128[pi][static_cast<std::size_t>(e_desc[i])];
        }

        u64 best_b = 0;
        auto dfs = [&](auto&& self, const int idx, const u64 cur) -> void {
            if (cur > limit) return;
            if (idx == pcount) {
                if (cur > best_b) best_b = cur;
                return;
            }
            if (static_cast<u128>(cur) * rem_max[idx] <= static_cast<u128>(best_b)) return;

            const int pi = desc_idx[idx];
            const auto& pw = pow_u64[pi];
            const int e = e_desc[idx];
            const u64 lim = limit / cur;
            int kmax = e;
            while (kmax > 0 && pw[static_cast<std::size_t>(kmax)] > lim) --kmax;

            for (int k = kmax; k >= 0; --k) {
                const u64 nxt = cur * pw[static_cast<std::size_t>(k)];
                if (static_cast<u128>(nxt) * rem_max[idx + 1] <= static_cast<u128>(best_b)) break;
                self(self, idx + 1, nxt);
            }
        };
        dfs(dfs, 0, 1ULL);
        return best_b;
    }

    void evaluate_candidate(const u64 a, const std::array<std::uint8_t, MAXP>& a_exp) {
        cpp_int qa = n_value / a;
        const u64 u = sqrt_floor_cpp(qa);

        std::array<int, MAXP> residual{};
        for (int i = 0; i < pcount; ++i) residual[i] = exps[i] - static_cast<int>(a_exp[i]);
        const u64 b = max_divisor_leq(residual, u);
        if (b < a) return;

        cpp_int c = qa / b;
        if (c < b) return;

        const cpp_int lhs = c * best.a;
        const cpp_int rhs = best.c * a;
        if (lhs < rhs) {
            best.a = a;
            best.b = b;
            best.c = c;
            best.sum = cpp_int(a) + cpp_int(b) + c;
            best_log_range = std::log(best.c.convert_to<long double>()) - std::log(best.a.convert_to<long double>());
        } else if (lhs == rhs) {
            const cpp_int s = cpp_int(a) + cpp_int(b) + c;
            if (s < best.sum) {
                best.a = a;
                best.b = b;
                best.c = c;
                best.sum = s;
            }
        }
    }

    void phase_b_search() {
        const u64 a_high = cube_root_floor();
        const u64 base_a_low = lower_a_bound(a_high);
        const long double base_log_low = std::log(static_cast<long double>(base_a_low));
        const long double log_high = std::log(static_cast<long double>(a_high));

        std::array<std::uint8_t, MAXP> a_exp{};
        a_exp.fill(0U);

        auto dfs = [&](auto&& self, const int idx, const u64 cur, const long double lv) -> void {
            const long double dynamic_log_low =
                std::max(base_log_low, (log_n - 2.0L * best_log_range) / 3.0L);
            if (idx == pcount) {
                if (lv + 1e-18L >= dynamic_log_low && cur <= a_high) {
                    evaluate_candidate(cur, a_exp);
                }
                return;
            }
            if (lv > log_high + 1e-18L) return;
            if (lv + rem_log[static_cast<std::size_t>(idx)] < dynamic_log_low - 1e-18L) return;

            const u64 p = primes[idx];
            const int e = exps[idx];
            u64 mul = 1;
            for (int k = 0; k <= e; ++k) {
                if (cur > a_high / mul) break;
                const u64 nxt = cur * mul;
                a_exp[idx] = static_cast<std::uint8_t>(k);
                self(self, idx + 1, nxt, lv + static_cast<long double>(k) * log_p[idx]);
                if (k < e) mul *= p;
            }
            a_exp[idx] = 0U;
        };
        dfs(dfs, 0, 1ULL, 0.0L);
    }

    cpp_int solve() {
        greedy_seed();

        std::array<long double, 3> logs = {0.0L, 0.0L, 0.0L};
        std::array<std::array<std::uint8_t, MAXP>, 3> bins{};
        for (auto& row : bins) row.fill(0U);
        phase_a_nodes = 0;
        phase_a_dfs(0, logs, bins);

        phase_b_search();
        return best.sum;
    }
};

std::pair<std::vector<u64>, std::vector<int>> factorize_number(u64 n) {
    std::vector<u64> primes;
    std::vector<int> exps;
    for (u64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) continue;
        int e = 0;
        while (n % p == 0) {
            n /= p;
            ++e;
        }
        primes.push_back(p);
        exps.push_back(e);
    }
    if (n > 1) {
        primes.push_back(n);
        exps.push_back(1);
    }
    return {primes, exps};
}

std::pair<std::vector<u64>, std::vector<int>> factorize_factorial(const int n) {
    std::vector<u64> primes;
    std::vector<int> exps;
    std::vector<bool> sieve(static_cast<std::size_t>(n + 1), true);
    sieve[0] = false;
    if (n >= 1) sieve[1] = false;
    for (int i = 2; i <= n; ++i) {
        if (!sieve[static_cast<std::size_t>(i)]) continue;
        for (int j = i + i; j <= n; j += i) sieve[static_cast<std::size_t>(j)] = false;
        primes.push_back(static_cast<u64>(i));
        int e = 0;
        int t = n;
        while (t > 0) {
            t /= i;
            e += t;
        }
        exps.push_back(e);
    }
    return {primes, exps};
}

bool run_checkpoints() {
    {
        auto [p, e] = factorize_number(165ULL);
        Solver solver(p, e);
        if (solver.solve() != cpp_int(19)) return false;
    }
    {
        auto [p, e] = factorize_number(100100ULL);
        Solver solver(p, e);
        if (solver.solve() != cpp_int(142)) return false;
    }
    {
        auto [p, e] = factorize_factorial(20);
        Solver solver(p, e);
        if (solver.solve() != cpp_int(4034872)) return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        std::cerr << "Checkpoint failed\n";
        return 2;
    }

    auto [primes, exps] = factorize_factorial(43);
    Solver solver(primes, exps);
    std::cout << to_string_cpp(solver.solve()) << '\n';
    return 0;
}
