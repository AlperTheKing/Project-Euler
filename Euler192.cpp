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
using u128 = unsigned __int128;
using i128 = __int128;

constexpr u64 kDenominatorBound = 1'000'000'000'000ULL;
constexpr u64 kMaxN = 100'000ULL;

u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) <= x / (r + 1)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

bool is_square(u64 n) {
    const u64 r = isqrt_u64(n);
    return r * r == n;
}

u64 best_denominator_cf(u64 n, u64 d_bound) {
    const u64 a0 = isqrt_u64(n);
    if (a0 * a0 == n) {
        return 0;
    }

    u64 m = 0;
    u64 den = 1;
    u64 a = a0;

    u128 p_prev2 = 0;
    u128 q_prev2 = 1;
    u128 p_prev1 = 1;
    u128 q_prev1 = 0;

    while (true) {
        const u128 p = static_cast<u128>(a) * p_prev1 + p_prev2;
        const u128 q = static_cast<u128>(a) * q_prev1 + q_prev2;

        if (q > d_bound) {
            u128 best_q = q_prev1;

            if (d_bound > q_prev1) {
                const u128 t = (static_cast<u128>(d_bound) - q_prev2) / q_prev1;
                if (t > 0) {
                    const u128 q_semiconv = t * q_prev1 + q_prev2;
                    const i128 rhs = static_cast<i128>(den) *
                                         (2 * static_cast<i128>(t) * static_cast<i128>(q_prev1) +
                                          static_cast<i128>(q_prev2)) -
                                     static_cast<i128>(q_prev1) * static_cast<i128>(m);

                    if (rhs > 0) {
                        const i128 upper =
                            static_cast<i128>(q_prev1) * static_cast<i128>(a0 + 1);

                        if (rhs > upper) {
                            best_q = q_semiconv;
                        } else {
                            const i128 lhs_sq = static_cast<i128>(n) *
                                                static_cast<i128>(q_prev1) *
                                                static_cast<i128>(q_prev1);
                            const i128 rhs_sq = rhs * rhs;
                            if (rhs_sq > lhs_sq) {
                                best_q = q_semiconv;
                            }
                        }
                    }
                }
            }

            return static_cast<u64>(best_q);
        }

        p_prev2 = p_prev1;
        q_prev2 = q_prev1;
        p_prev1 = p;
        q_prev1 = q;

        const u64 m_next = den * a - m;
        const u64 den_next = (n - m_next * m_next) / den;
        const u64 a_next = (a0 + m_next) / den_next;

        m = m_next;
        den = den_next;
        a = a_next;
    }
}

u64 best_denominator_bruteforce(u64 n, u64 d_bound) {
    const long double sq = std::sqrt(static_cast<long double>(n));

    long double best_err = 0.0L;
    u64 best_den = 0;
    bool found = false;

    for (u64 q = 1; q <= d_bound; ++q) {
        const long double v = sq * static_cast<long double>(q);
        const u64 p0 = static_cast<u64>(std::floor(v));

        for (u64 p : {p0, p0 + 1}) {
            if (std::gcd(p, q) != 1) {
                continue;
            }

            const long double err = std::fabsl(static_cast<long double>(p) / static_cast<long double>(q) - sq);
            if (!found || err < best_err) {
                found = true;
                best_err = err;
                best_den = q;
            }
        }
    }

    return best_den;
}

bool run_validations() {
    if (best_denominator_cf(13, 20) != 5) {
        std::cerr << "Validation failed: sqrt(13), d=20 should yield denominator 5.\n";
        return false;
    }

    if (best_denominator_cf(13, 30) != 28) {
        std::cerr << "Validation failed: sqrt(13), d=30 should yield denominator 28.\n";
        return false;
    }

    for (u64 n = 2; n <= 200; ++n) {
        if (is_square(n)) {
            continue;
        }

        for (u64 d = 2; d <= 200; ++d) {
            const u64 fast = best_denominator_cf(n, d);
            const u64 brute = best_denominator_bruteforce(n, d);
            if (fast != brute) {
                std::cerr << "Validation failed: n=" << n << ", d=" << d
                          << ", fast=" << fast << ", brute=" << brute << "\n";
                return false;
            }
        }
    }

    return true;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string s;
    while (value > 0) {
        const u64 digit = static_cast<u64>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

u128 solve_parallel(u64 max_n, u64 d_bound, unsigned thread_count) {
    if (thread_count == 0) {
        thread_count = 1;
    }

    std::atomic<u64> next_n{2};
    std::vector<u128> partial(thread_count, 0);
    std::vector<std::thread> workers;
    workers.reserve(thread_count);

    for (unsigned tid = 0; tid < thread_count; ++tid) {
        workers.emplace_back([&, tid]() {
            u128 local_sum = 0;

            while (true) {
                const u64 n = next_n.fetch_add(1, std::memory_order_relaxed);
                if (n > max_n) {
                    break;
                }

                if (is_square(n)) {
                    continue;
                }

                local_sum += best_denominator_cf(n, d_bound);
            }

            partial[tid] = local_sum;
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }

    u128 total = 0;
    for (u128 part : partial) {
        total += part;
    }
    return total;
}

}  // namespace

int main() {
    if (!run_validations()) {
        return 1;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }

    const u128 answer = solve_parallel(kMaxN, kDenominatorBound, threads);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
