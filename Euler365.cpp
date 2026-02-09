#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

i64 mod_pow(i64 base, i64 exp, i64 mod) {
    i64 result = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if ((exp & 1LL) != 0LL) {
            result = static_cast<i64>((static_cast<__int128>(result) * base) % mod);
        }
        base = static_cast<i64>((static_cast<__int128>(base) * base) % mod);
        exp >>= 1LL;
    }
    return result;
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }
    std::string out;
    while (value > 0U) {
        const unsigned digit = static_cast<unsigned>(value % 10U);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

std::vector<int> primes_between(const int lo_exclusive, const int hi_exclusive) {
    std::vector<bool> is_prime(static_cast<std::size_t>(hi_exclusive), true);
    if (hi_exclusive > 0) {
        is_prime[0] = false;
    }
    if (hi_exclusive > 1) {
        is_prime[1] = false;
    }
    for (int p = 2; p * p < hi_exclusive; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q < hi_exclusive; q += p) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }

    std::vector<int> primes;
    for (int x = std::max(2, lo_exclusive + 1); x < hi_exclusive; ++x) {
        if (is_prime[static_cast<std::size_t>(x)]) {
            primes.push_back(x);
        }
    }
    return primes;
}

int binom_mod_prime_lucas(const u64 n, const u64 k, const int p) {
    std::vector<int> fact(static_cast<std::size_t>(p), 1);
    for (int i = 1; i < p; ++i) {
        fact[static_cast<std::size_t>(i)] =
            static_cast<int>((static_cast<i64>(fact[static_cast<std::size_t>(i - 1)]) * i) % p);
    }

    std::vector<int> inv_fact(static_cast<std::size_t>(p), 1);
    inv_fact[static_cast<std::size_t>(p - 1)] =
        static_cast<int>(mod_pow(fact[static_cast<std::size_t>(p - 1)], p - 2, p));
    for (int i = p - 1; i >= 1; --i) {
        inv_fact[static_cast<std::size_t>(i - 1)] =
            static_cast<int>((static_cast<i64>(inv_fact[static_cast<std::size_t>(i)]) * i) % p);
    }

    u64 nn = n;
    u64 kk = k;
    i64 result = 1;
    while (nn > 0 || kk > 0) {
        const int ni = static_cast<int>(nn % static_cast<u64>(p));
        const int ki = static_cast<int>(kk % static_cast<u64>(p));
        if (ki > ni) {
            return 0;
        }

        i64 term = fact[static_cast<std::size_t>(ni)];
        term = (term * inv_fact[static_cast<std::size_t>(ki)]) % p;
        term = (term * inv_fact[static_cast<std::size_t>(ni - ki)]) % p;
        result = (result * term) % p;

        nn /= static_cast<u64>(p);
        kk /= static_cast<u64>(p);
    }

    return static_cast<int>(result);
}

u128 sum_crt_binom_over_triples(const u64 n, const u64 k, const std::vector<int>& primes) {
    const int m = static_cast<int>(primes.size());
    std::vector<int> residues(static_cast<std::size_t>(m), 0);
    for (int i = 0; i < m; ++i) {
        residues[static_cast<std::size_t>(i)] = binom_mod_prime_lucas(n, k, primes[static_cast<std::size_t>(i)]);
    }

    std::vector<std::vector<int>> inv(static_cast<std::size_t>(m), std::vector<int>(static_cast<std::size_t>(m), 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i == j) {
                continue;
            }
            const int mod = primes[static_cast<std::size_t>(j)];
            inv[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                static_cast<int>(mod_pow(primes[static_cast<std::size_t>(i)] % mod, mod - 2, mod));
        }
    }

    u128 total = 0;
    for (int i = 0; i < m - 2; ++i) {
        const i64 p = primes[static_cast<std::size_t>(i)];
        const i64 a = residues[static_cast<std::size_t>(i)];

        for (int j = i + 1; j < m - 1; ++j) {
            const i64 q = primes[static_cast<std::size_t>(j)];
            const i64 b = residues[static_cast<std::size_t>(j)];

            const i64 diff_ab = (b - a + q) % q;
            const i64 t = (diff_ab * inv[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]) % q;
            const i64 x_pq = a + p * t;  // modulo p*q
            const i64 pq = p * q;

            for (int kidx = j + 1; kidx < m; ++kidx) {
                const i64 r = primes[static_cast<std::size_t>(kidx)];
                const i64 c = residues[static_cast<std::size_t>(kidx)];

                const i64 diff_c = (c - (x_pq % r) + r) % r;
                const i64 inv_pq_mod_r =
                    (static_cast<i64>(inv[static_cast<std::size_t>(i)][static_cast<std::size_t>(kidx)]) *
                     inv[static_cast<std::size_t>(j)][static_cast<std::size_t>(kidx)]) %
                    r;
                const i64 u = (diff_c * inv_pq_mod_r) % r;
                const i64 x = x_pq + pq * u;  // modulo p*q*r

                total += static_cast<u128>(x);
            }
        }
    }

    return total;
}

u128 binom_exact_small(const int n, const int k) {
    const int kk = std::min(k, n - k);
    std::vector<u64> num;
    num.reserve(static_cast<std::size_t>(kk));
    for (int x = n - kk + 1; x <= n; ++x) {
        num.push_back(static_cast<u64>(x));
    }

    for (int d = 2; d <= kk; ++d) {
        u64 rem = static_cast<u64>(d);
        for (u64& v : num) {
            const u64 g = std::gcd(v, rem);
            if (g > 1U) {
                v /= g;
                rem /= g;
                if (rem == 1U) {
                    break;
                }
            }
        }
    }

    u128 result = 1;
    for (const u64 v : num) {
        result *= static_cast<u128>(v);
    }
    return result;
}

bool run_checkpoints() {
    // Lucas check against small direct binomial modulo prime.
    const u128 exact_30_12 = binom_exact_small(30, 12);
    for (const int p : {7, 11, 13, 17, 19, 23, 29}) {
        const int lucas = binom_mod_prime_lucas(30ULL, 12ULL, p);
        const int direct = static_cast<int>(exact_30_12 % static_cast<u128>(p));
        if (lucas != direct) {
            std::cerr << "Checkpoint failed: Lucas mismatch for p=" << p << '\n';
            return false;
        }
    }

    // CRT+triple sum check on a tiny instance with exact arithmetic.
    const std::vector<int> tiny_primes = {7, 11, 13, 17};
    const u128 exact_20_8 = binom_exact_small(20, 8);
    u128 direct_sum = 0;
    for (int i = 0; i < static_cast<int>(tiny_primes.size()) - 2; ++i) {
        for (int j = i + 1; j < static_cast<int>(tiny_primes.size()) - 1; ++j) {
            for (int k = j + 1; k < static_cast<int>(tiny_primes.size()); ++k) {
                const u64 mod = static_cast<u64>(tiny_primes[static_cast<std::size_t>(i)]) *
                                static_cast<u64>(tiny_primes[static_cast<std::size_t>(j)]) *
                                static_cast<u64>(tiny_primes[static_cast<std::size_t>(k)]);
                direct_sum += (exact_20_8 % static_cast<u128>(mod));
            }
        }
    }
    const u128 crt_sum = sum_crt_binom_over_triples(20ULL, 8ULL, tiny_primes);
    if (crt_sum != direct_sum) {
        std::cerr << "Checkpoint failed: CRT triple sum mismatch\n";
        return false;
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
        return 2;
    }

    const std::vector<int> primes = primes_between(1000, 5000);
    const u128 answer = sum_crt_binom_over_triples(1000000000000000000ULL, 1000000000ULL, primes);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
