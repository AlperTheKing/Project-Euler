#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;

class Solver {
public:
    explicit Solver(u64 n) : N_(n) {
        rN_ = static_cast<int>(std::sqrt(static_cast<long double>(N_)));
        while (static_cast<u64>(rN_ + 1) * static_cast<u64>(rN_ + 1) <= N_) {
            ++rN_;
        }
        while (static_cast<u64>(rN_) * static_cast<u64>(rN_) > N_) {
            --rN_;
        }
        primes_.push_back(2);
        get_primes(2 * rN_);
        pc_ = lucy();
    }

    u64 Q(int k) const {
        if (k < 2) {
            return 0;
        }
        return cumulative(k) - cumulative(k - 1);
    }

    u64 cumulative(int k) const {
        if (k < 2) {
            return 0;
        }
        if (k < static_cast<int>(cumulative_cache_.size()) &&
            cumulative_cache_[static_cast<std::size_t>(k)] != std::numeric_limits<u64>::max()) {
            return cumulative_cache_[static_cast<std::size_t>(k)];
        }

        const int two_k = 2 * k;
        u64 ans = 0;
        std::vector<int> cp;

        auto power_of_two_index = [](u64 x) -> int {
            if (x == 1) return 0;
            if (x == 2) return 1;
            if (x == 4) return 2;
            if (x == 8) return 3;
            if (x == 16) return 4;
            return -1;
        };

        dfs(0, N_, 1, 1, cp, two_k, [&](u64 pw, u64 cb, const std::vector<int>& cur_cp) {
            for (int k0 = 2; k0 <= k; ++k0) {
                const u64 two_k0 = static_cast<u64>(2 * k0);
                if (two_k0 % cb != 0) {
                    continue;
                }
                const int k2 = power_of_two_index(two_k0 / cb);
                if (k2 < 0) {
                    continue;
                }
                if (k2 == 0) {
                    ans += 1;
                } else {
                    ans += rec(cur_cp, 0, N_ / pw, pw, k2);
                }
            }
        });

        if (k < static_cast<int>(cumulative_cache_.size())) {
            cumulative_cache_[static_cast<std::size_t>(k)] = ans;
        }
        return ans;
    }

private:
    u64 N_;
    int rN_;
    std::vector<int> primes_;
    std::vector<i64> pc_;
    mutable std::vector<u64> cumulative_cache_{32, std::numeric_limits<u64>::max()};

    static bool contains_prime(const std::vector<int>& cp, int p) {
        for (int q : cp) {
            if (q == p) {
                return true;
            }
        }
        return false;
    }

    i64 pc_at(int idx) const {
        if (idx >= 0) {
            return pc_[static_cast<std::size_t>(idx)];
        }
        return pc_[pc_.size() - static_cast<std::size_t>(-idx)];
    }

    void get_primes(int limit) {
        std::vector<int> sieve(limit + 1, 0);
        for (int i = 2; i <= limit; i += 2) {
            sieve[static_cast<std::size_t>(i)] = 2;
        }

        for (int i = 3; i <= limit; i += 2) {
            if (sieve[static_cast<std::size_t>(i)] == 0) {
                primes_.push_back(i);
                sieve[static_cast<std::size_t>(i)] = i;
                if (static_cast<i64>(i) * static_cast<i64>(i) <= limit) {
                    for (int j = i * i; j <= limit; j += 2 * i) {
                        if (sieve[static_cast<std::size_t>(j)] == 0) {
                            sieve[static_cast<std::size_t>(j)] = i;
                        }
                    }
                }
            }
        }
    }

    std::vector<i64> lucy() const {
        const int r = rN_;
        std::vector<i64> S;
        S.reserve(static_cast<std::size_t>(2 * r + 1));
        for (int i = 0; i <= r; ++i) {
            S.push_back(static_cast<i64>(i - 1));
        }
        for (int i = r; i >= 1; --i) {
            S.push_back(static_cast<i64>(N_ / static_cast<u64>(i) - 1));
        }

        auto at = [&](int idx) -> i64& {
            if (idx >= 0) {
                return S[static_cast<std::size_t>(idx)];
            }
            return S[S.size() - static_cast<std::size_t>(-idx)];
        };

        for (int p = 2; p <= r; ++p) {
            if (S[static_cast<std::size_t>(p)] <= S[static_cast<std::size_t>(p - 1)]) {
                continue;
            }

            const i64 sp = S[static_cast<std::size_t>(p - 1)];
            const u64 p2 = static_cast<u64>(p) * static_cast<u64>(p);

            for (int i = 1; i <= r; ++i) {
                if (N_ / static_cast<u64>(i) < p2) {
                    break;
                }
                const u64 ip = static_cast<u64>(i) * static_cast<u64>(p);
                const u64 nip = N_ / ip;
                const int idx = (nip <= static_cast<u64>(r)) ? static_cast<int>(nip) : -static_cast<int>(ip);
                at(-i) -= (at(idx) - sp);
            }

            for (int i = r; i >= 1; --i) {
                if (static_cast<u64>(i) < p2) {
                    break;
                }
                S[static_cast<std::size_t>(i)] -= (S[static_cast<std::size_t>(i / p)] - sp);
            }
        }

        return S;
    }

    u64 rec(const std::vector<int>& cp, int i, u64 n, u64 cn, int k) const {
        if (k == 1) {
            const int idx = (cn < static_cast<u64>(rN_)) ? -static_cast<int>(cn) : static_cast<int>(N_ / cn);
            i64 ans = pc_at(idx) - i;
            if (ans < 0) {
                ans = 0;
            }
            const int start_prime =
                (i < static_cast<int>(primes_.size())) ? primes_[static_cast<std::size_t>(i)] : std::numeric_limits<int>::max();
            for (int p : cp) {
                if (static_cast<u64>(p) <= n && p >= start_prime) {
                    --ans;
                }
            }
            return (ans > 0) ? static_cast<u64>(ans) : 0;
        }

        u64 ans = 0;
        while (true) {
            if (i >= static_cast<int>(primes_.size())) {
                break;
            }
            const int p = primes_[static_cast<std::size_t>(i)];
            if (contains_prime(cp, p)) {
                ++i;
                continue;
            }
            const u64 nn = n / static_cast<u64>(p);
            if (nn < static_cast<u64>(p)) {
                break;
            }
            ans += rec(cp, i + 1, nn, cn * static_cast<u64>(p), k - 1);
            ++i;
        }
        return ans;
    }

    template <class F>
    void dfs(int i, u64 n, u64 c, u64 cb, std::vector<int>& cp, int two_k, F&& emit) const {
        emit(c, cb, cp);

        while (i < static_cast<int>(primes_.size())) {
            const u64 p = static_cast<u64>(primes_[static_cast<std::size_t>(i)]);
            if (p > std::numeric_limits<u64>::max() / p) {
                break;
            }

            u64 nn = n / (p * p);
            if (nn == 0) {
                break;
            }

            u64 cc = c * p * p;
            int e = 2;
            cp.push_back(static_cast<int>(p));

            while (true) {
                const u64 mult = static_cast<u64>(e + (e & 1));
                if (cb <= static_cast<u64>(two_k) / mult) {
                    dfs(i + 1, nn, cc, cb * mult, cp, two_k, emit);
                }

                if (nn < p) {
                    break;
                }
                nn /= p;
                if (cc > std::numeric_limits<u64>::max() / p) {
                    break;
                }
                cc *= p;
                ++e;
            }

            cp.pop_back();
            ++i;
        }
    }
};

int main() {
    {
        Solver s(100);
        assert(s.Q(2) == 51);
    }
    {
        Solver s(1'000'000);
        assert(s.Q(6) == 6'189);
    }

    Solver s(1'000'000'000'000ULL);
    u64 ans = 0;
    for (int k = 2; k <= 10; ++k) {
        ans += s.Q(k);
    }
    std::cout << ans << '\n';
    return 0;
}
