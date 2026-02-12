#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <iostream>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr i64 kMod = 1'234'567'891;

struct SimResult {
    i64 end_t;
    std::vector<std::vector<int>> patterns;
    int kmax;
};

static std::vector<int> sieve_spf(int n) {
    std::vector<int> spf(n + 1);
    for (int i = 0; i <= n; ++i) {
        spf[i] = i;
    }
    int lim = static_cast<int>(std::sqrt(static_cast<long double>(n)));
    for (int i = 2; i <= lim; ++i) {
        if (spf[i] == i) {
            for (int j = i * i; j <= n; j += i) {
                if (spf[j] == j) {
                    spf[j] = i;
                }
            }
        }
    }
    return spf;
}

static std::vector<int> factor_list(int x, const std::vector<int>& spf) {
    std::vector<int> out;
    while (x > 1) {
        int p = spf[x];
        out.push_back(p);
        x /= p;
    }
    return out;
}

static i64 direct_sum_mod(int n, i64 m, i64 mod) {
    const auto spf = sieve_spf(n);

    std::vector<std::vector<int>> piles;
    std::vector<int> pos;
    piles.reserve(n);
    pos.reserve(n);
    for (int i = 2; i <= n; ++i) {
        piles.push_back(factor_list(i, spf));
        pos.push_back(0);
    }

    for (i64 t = 0; t < m; ++t) {
        const int k = static_cast<int>(piles.size());
        std::vector<int> extracted(k);
        std::vector<std::vector<int>> next_piles;
        std::vector<int> next_pos;
        next_piles.reserve(k + 1);
        next_pos.reserve(k + 1);

        for (int i = 0; i < k; ++i) {
            extracted[i] = piles[i][pos[i]];
            ++pos[i];
            if (pos[i] < static_cast<int>(piles[i].size())) {
                next_piles.push_back(std::move(piles[i]));
                next_pos.push_back(pos[i]);
            }
        }

        std::sort(extracted.begin(), extracted.end());
        next_piles.push_back(std::move(extracted));
        next_pos.push_back(0);

        piles.swap(next_piles);
        pos.swap(next_pos);
    }

    i64 total = 0;
    for (int i = 0; i < static_cast<int>(piles.size()); ++i) {
        i64 prod = 1;
        for (int j = pos[i]; j < static_cast<int>(piles[i].size()); ++j) {
            prod = static_cast<i64>((static_cast<u128>(prod) * static_cast<u64>(piles[i][j])) % mod);
        }
        total += prod;
        if (total >= mod) {
            total -= mod;
        }
    }
    return total;
}

static SimResult simulate_until_periodic(
    int n,
    int k_extra = 10,
    int streak_needed = 2000,
    int max_rounds = 200000
) {
    const auto spf = sieve_spf(n);

    std::vector<std::vector<int>> piles;
    std::vector<int> pos;
    piles.reserve(n);
    pos.reserve(n);

    i64 total_factors = 0;
    for (int i = 2; i <= n; ++i) {
        auto f = factor_list(i, spf);
        total_factors += static_cast<i64>(f.size());
        piles.push_back(std::move(f));
        pos.push_back(0);
    }

    const int k_lim = static_cast<int>(std::sqrt(2.0L * static_cast<long double>(total_factors))) + k_extra;

    std::vector<std::deque<int>> bufs(k_lim + 1);

    int stable_streak = 0;

    for (int t = 1; t <= max_rounds; ++t) {
        const int k = static_cast<int>(piles.size());
        std::vector<int> extracted(k);

        std::vector<std::vector<int>> next_piles;
        std::vector<int> next_pos;
        next_piles.reserve(k + 1);
        next_pos.reserve(k + 1);

        for (int i = 0; i < k; ++i) {
            extracted[i] = piles[i][pos[i]];
            ++pos[i];
            if (pos[i] < static_cast<int>(piles[i].size())) {
                next_piles.push_back(std::move(piles[i]));
                next_pos.push_back(pos[i]);
            }
        }

        std::sort(extracted.begin(), extracted.end());
        next_piles.push_back(extracted);
        next_pos.push_back(0);

        piles.swap(next_piles);
        pos.swap(next_pos);

        int mm = std::min(static_cast<int>(extracted.size()), k_lim);

        if (t <= k_lim) {
            for (int kk = 1; kk <= mm; ++kk) {
                auto& b = bufs[kk];
                b.push_back(extracted[kk - 1]);
                if (static_cast<int>(b.size()) > kk) {
                    b.pop_front();
                }
            }
            for (int kk = mm + 1; kk <= k_lim; ++kk) {
                auto& b = bufs[kk];
                b.push_back(1);
                if (static_cast<int>(b.size()) > kk) {
                    b.pop_front();
                }
            }
            stable_streak = 0;
            continue;
        }

        bool all_ok = true;

        for (int kk = 1; kk <= mm; ++kk) {
            auto& b = bufs[kk];
            const int v = extracted[kk - 1];
            if (b.front() != v) {
                all_ok = false;
            }
            b.push_back(v);
            if (static_cast<int>(b.size()) > kk) {
                b.pop_front();
            }
        }

        for (int kk = mm + 1; kk <= k_lim; ++kk) {
            auto& b = bufs[kk];
            if (b.front() != 1) {
                all_ok = false;
            }
            b.push_back(1);
            if (static_cast<int>(b.size()) > kk) {
                b.pop_front();
            }
        }

        if (all_ok) {
            ++stable_streak;
            if (stable_streak >= streak_needed) {
                std::vector<std::vector<int>> patterns(k_lim + 1);
                int kmax = 0;
                for (int kk = 1; kk <= k_lim; ++kk) {
                    patterns[kk] = std::vector<int>(bufs[kk].begin(), bufs[kk].end());
                    bool has_non_one = false;
                    for (int x : patterns[kk]) {
                        if (x != 1) {
                            has_non_one = true;
                            break;
                        }
                    }
                    if (has_non_one) {
                        kmax = kk;
                    }
                }
                return {t, std::move(patterns), kmax};
            }
        } else {
            stable_streak = 0;
        }
    }

    throw std::runtime_error("periodicity not detected");
}

static i64 sum_at_round_mod(int n, i64 m, i64 mod) {
    const SimResult sim = simulate_until_periodic(n);

    if (m <= sim.end_t) {
        return direct_sum_mod(n, m, mod);
    }

    const i64 r0 = m - sim.end_t - 1;
    i64 total = 0;

    for (int d = 0; d < sim.kmax; ++d) {
        const i64 r = r0 - d;
        if (sim.patterns[d + 1][static_cast<int>(r % (d + 1))] == 1) {
            continue;
        }

        i64 prod = 1;
        for (int k = sim.kmax; k > d; --k) {
            const int v = sim.patterns[k][static_cast<int>(r % k)];
            if (v != 1) {
                prod = static_cast<i64>((static_cast<u128>(prod) * static_cast<u64>(v)) % mod);
            }
        }

        total += prod;
        if (total >= mod) {
            total -= mod;
        }
    }

    return total;
}

int main() {
    if (direct_sum_mod(5, 3, (1LL << 62)) != 21) {
        std::cerr << "Validation failed for S(5,3).\n";
        return 1;
    }
    if (direct_sum_mod(10, 100, (1LL << 62)) != 257) {
        std::cerr << "Validation failed for S(10,100).\n";
        return 1;
    }

    std::cout << sum_at_round_mod(10'000, 10'000'000'000'000'000LL, kMod) << '\n';
    return 0;
}
