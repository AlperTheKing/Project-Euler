#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;

u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) <= x / (r + 1)) ++r;
    while (r > x / r) --r;
    return r;
}

class Solver578 {
public:
    explicit Solver578(u64 n) : n_(n), r_(isqrt_u64(n)), b_(n / r_) {
        std::vector<u64> values;
        values.reserve(static_cast<std::size_t>(r_ + b_));
        for (u64 v = 1; v <= r_; ++v) values.push_back(n_ / v);
        for (u64 v = b_; v >= 1; --v) {
            if (v == b_) continue;
            values.push_back(v);
        }

        count0_.resize(static_cast<std::size_t>(b_));
        for (u64 v = 1; v <= b_; ++v) count0_[static_cast<std::size_t>(v - 1)] = v - 1;

        count1_.resize(static_cast<std::size_t>(r_));
        for (u64 v = 1; v <= r_; ++v) count1_[static_cast<std::size_t>(v - 1)] = n_ / v - 1;

        for (u64 p = 2; p <= r_; ++p) {
            const u64 prev = count0_[static_cast<std::size_t>(p - 2)];
            if (count0_[static_cast<std::size_t>(p - 1)] <= prev) continue;

            primes_.push_back(p);
            const u64 p2 = p * p;

            for (u64 v : values) {
                if (v < p2) break;

                if (v <= b_) {
                    count0_[static_cast<std::size_t>(v - 1)] -=
                        count0_[static_cast<std::size_t>(v / p - 1)] - prev;
                } else {
                    const std::size_t idx = static_cast<std::size_t>(n_ / v - 1);
                    const u64 vp = v / p;
                    if (vp <= b_) {
                        count1_[idx] -= count0_[static_cast<std::size_t>(vp - 1)] - prev;
                    } else {
                        count1_[idx] -= count1_[static_cast<std::size_t>(n_ / vp - 1)] - prev;
                    }
                }
            }
        }

        prime_pows_.resize(primes_.size());
        for (std::size_t idx = 0; idx < primes_.size(); ++idx) {
            const u64 p = primes_[idx];
            auto& pw = prime_pows_[idx];
            pw.reserve(8);
            pw.push_back(1);
            while (true) {
                const u64 cur = pw.back();
                if (cur > n_ / p) break;
                pw.push_back(cur * p);
            }
        }
    }

    u64 solve(unsigned int thread_count = 0) const {
        if (primes_.empty()) return 1;
        std::vector<std::pair<unsigned char, u64>> roots;
        roots.reserve(64);
        const u64 first_prime = primes_[0];
        u64 p_pow = first_prime;
        for (int e = 1; e <= 100 && p_pow <= n_; ++e) {
            roots.emplace_back(static_cast<unsigned char>(e), p_pow);
            if (p_pow > n_ / first_prime) break;
            p_pow *= first_prime;
        }

        if (roots.empty()) return 1;

        if (thread_count == 0) {
            thread_count = std::thread::hardware_concurrency();
            if (thread_count == 0) thread_count = 1;
        }
        thread_count = std::min<unsigned int>(thread_count, static_cast<unsigned int>(roots.size()));

        if (thread_count <= 1) {
            return 1 + solve_roots_sequential(roots);
        }

        std::atomic<std::size_t> next(0);
        std::vector<u64> partial(static_cast<std::size_t>(thread_count), 0);
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(thread_count));

        for (unsigned int t = 0; t < thread_count; ++t) {
            workers.emplace_back([&, t]() {
                std::vector<unsigned char> current;
                current.reserve(64);
                u64 local = 0;
                while (true) {
                    const std::size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                    if (idx >= roots.size()) break;
                    const auto [e, pw] = roots[idx];
                    current.clear();
                    current.push_back(e);
                    local += signature_count(current, 0, n_, 0);
                    local += enumerate_sum(n_ / pw, 1, e, current);
                }
                partial[static_cast<std::size_t>(t)] = local;
            });
        }
        for (auto& th : workers) th.join();
        u64 total = 1;
        for (u64 v : partial) total += v;
        return total;
    }

private:
    u64 n_;
    u64 r_;
    u64 b_;
    std::vector<u64> primes_;
    std::vector<u64> count0_;
    std::vector<u64> count1_;
    std::vector<std::vector<u64>> prime_pows_;

    bool pow_leq_idx(std::size_t idx, int e, u64 limit, u64& out) const {
        if (e < 0) return false;
        const auto& pw = prime_pows_[idx];
        if (static_cast<std::size_t>(e) >= pw.size()) return false;
        const u64 val = pw[static_cast<std::size_t>(e)];
        if (val > limit) return false;
        out = val;
        return true;
    }

    u64 pi(u64 v) const {
        if (v <= 1) return 0;
        if (v <= b_) return count0_[static_cast<std::size_t>(v - 1)];
        return count1_[static_cast<std::size_t>(n_ / v - 1)];
    }

    u64 count_pairs(int e0, int e1, u64 limit, std::size_t i) const {
        if (i + 1 >= primes_.size()) return 0;

        u64 total = 0;
        for (std::size_t j = i; j < primes_.size(); ++j) {
            u64 p0e = 1;
            if (!pow_leq_idx(j, e0, limit, p0e)) break;

            const u64 lim1 = limit / p0e;
            if (e1 == 1) {
                const u64 cnt = pi(lim1);
                if (cnt <= j + 1) break;
                total += cnt - (j + 1);
            } else if (e1 == 2) {
                const u64 cnt = pi(isqrt_u64(lim1));
                if (cnt <= j + 1) break;
                total += cnt - (j + 1);
            } else if (e1 == 4) {
                const u64 cnt = pi(isqrt_u64(isqrt_u64(lim1)));
                if (cnt <= j + 1) break;
                total += cnt - (j + 1);
            } else if (e1 == 8) {
                const u64 cnt = pi(isqrt_u64(isqrt_u64(isqrt_u64(lim1))));
                if (cnt <= j + 1) break;
                total += cnt - (j + 1);
            } else {
                for (std::size_t k = j + 1; k < primes_.size(); ++k) {
                    u64 p1e = 1;
                    if (!pow_leq_idx(k, e1, lim1, p1e)) break;
                    ++total;
                }
            }
        }

        return total;
    }

    u64 signature_count(const std::vector<unsigned char>& s, int pos, u64 limit,
                        std::size_t i) const {
        const int rem = static_cast<int>(s.size()) - pos;
        if (rem == 1) {
            const int e = static_cast<int>(s[pos]);
            if (e == 1) {
                const u64 cnt = pi(limit);
                return (cnt > i) ? cnt - i : 0ULL;
            }
            if (e == 2) {
                const u64 cnt = pi(isqrt_u64(limit));
                return (cnt > i) ? cnt - i : 0ULL;
            }
            if (e == 4) {
                const u64 cnt = pi(isqrt_u64(isqrt_u64(limit)));
                return (cnt > i) ? cnt - i : 0ULL;
            }
            if (e == 8) {
                const u64 cnt = pi(isqrt_u64(isqrt_u64(isqrt_u64(limit))));
                return (cnt > i) ? cnt - i : 0ULL;
            }

            u64 cnt = 0;
            for (std::size_t j = i; j < primes_.size(); ++j) {
                u64 pw = 1;
                if (!pow_leq_idx(j, e, limit, pw)) break;
                ++cnt;
            }
            return cnt;
        }
        if (rem == 2) {
            return count_pairs(static_cast<int>(s[pos]), static_cast<int>(s[pos + 1]), limit, i);
        }

        u64 total = 0;
        for (std::size_t j = i; j < primes_.size(); ++j) {
            bool ok = true;
            u64 minimal = 1;
            const std::size_t max_t =
                std::min<std::size_t>(static_cast<std::size_t>(rem), primes_.size() - j);

            for (std::size_t t = 0; t < max_t; ++t) {
                u64 pw = 1;
                if (!pow_leq_idx(j + t, static_cast<int>(s[pos + static_cast<int>(t)]),
                                 limit / minimal, pw)) {
                    ok = false;
                    break;
                }
                minimal *= pw;
            }
            if (!ok) break;

            u64 first_pw = 1;
            if (!pow_leq_idx(j, static_cast<int>(s[pos]), limit, first_pw)) break;

            total += signature_count(s, pos + 1, limit / first_pw, j + 1);
        }

        return total;
    }

    u64 solve_roots_sequential(const std::vector<std::pair<unsigned char, u64>>& roots) const {
        std::vector<unsigned char> current;
        current.reserve(64);
        u64 total = 0;
        for (const auto [e, pw] : roots) {
            current.clear();
            current.push_back(e);
            total += signature_count(current, 0, n_, 0);
            total += enumerate_sum(n_ / pw, 1, e, current);
        }
        return total;
    }

    u64 enumerate_sum(u64 limit, std::size_t i, int exponent_bound,
                      std::vector<unsigned char>& current) const {
        if (i >= primes_.size()) return 0;
        u64 total = 0;

        const u64 p = primes_[i];
        u64 p_pow = p;

        for (int e = 1; e <= exponent_bound && p_pow <= limit; ++e) {
            current.push_back(static_cast<unsigned char>(e));
            total += signature_count(current, 0, n_, 0);
            total += enumerate_sum(limit / p_pow, i + 1, e, current);
            current.pop_back();

            if (p_pow > limit / p) break;
            p_pow *= p;
        }
        return total;
    }
};

bool is_decreasing_prime_power(u64 n) {
    if (n <= 1) return true;
    std::vector<int> exponents;
    u64 x = n;

    for (u64 p = 2; p * p <= x; ++p) {
        if (x % p != 0) continue;
        int cnt = 0;
        while (x % p == 0) {
            x /= p;
            ++cnt;
        }
        exponents.push_back(cnt);
    }

    if (x > 1) exponents.push_back(1);

    for (std::size_t i = 1; i < exponents.size(); ++i) {
        if (exponents[i - 1] < exponents[i]) return false;
    }
    return true;
}

u64 brute_count(u64 n) {
    u64 cnt = 0;
    for (u64 i = 1; i <= n; ++i) {
        if (is_decreasing_prime_power(i)) ++cnt;
    }
    return cnt;
}

bool run_validations() {
    struct Check {
        u64 n;
        u64 expected;
    };

    const std::vector<Check> checks = {
        {100ULL, 94ULL},
        {1'000'000ULL, 922'052ULL},
    };

    for (const auto& chk : checks) {
        Solver578 solver(chk.n);
        const u64 got = solver.solve();
        if (got != chk.expected) {
            std::cerr << "Validation failed for C(" << chk.n << "): got " << got
                      << ", expected " << chk.expected << "\n";
            return false;
        }
    }

    const u64 small_n = 10'000ULL;
    Solver578 small_solver(small_n);
    const u64 got_small = small_solver.solve();
    const u64 expected_small = brute_count(small_n);
    if (got_small != expected_small) {
        std::cerr << "Validation failed for brute-force C(" << small_n << "): got "
                  << got_small << ", expected " << expected_small << "\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    u64 n = 10'000'000'000'000ULL;
    if (argc > 1) n = std::strtoull(argv[1], nullptr, 10);

    if (!run_validations()) return 1;

    Solver578 solver(n);
    std::cout << solver.solve() << '\n';
    return 0;
}
