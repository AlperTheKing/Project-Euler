#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int kLimit = 20000;
constexpr int kPrimeSearchLimit = 2000000;

std::vector<int> sieve_primes(int n) {
    std::vector<bool> is_prime(n + 1, true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (!is_prime[i]) {
            continue;
        }
        for (long long j = 1LL * i * i; j <= n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = false;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

std::vector<std::pair<int, int>> factorize(int n, const std::vector<int>& primes) {
    std::vector<std::pair<int, int>> factors;
    int x = n;
    for (int q : primes) {
        if (1LL * q * q > x) {
            break;
        }
        if (x % q != 0) {
            continue;
        }
        int e = 0;
        while (x % q == 0) {
            x /= q;
            ++e;
        }
        factors.push_back({q, e});
    }
    if (x > 1) {
        factors.push_back({x, 1});
    }
    return factors;
}

std::vector<int> divisors_from_factors(const std::vector<std::pair<int, int>>& factors) {
    std::vector<int> divisors{1};
    for (auto [prime, exponent] : factors) {
        const std::vector<int> previous = divisors;
        int power = 1;
        for (int e = 1; e <= exponent; ++e) {
            power *= prime;
            for (int d : previous) {
                divisors.push_back(d * power);
            }
        }
    }
    std::sort(divisors.begin(), divisors.end());
    return divisors;
}

int primitive_root(int p, const std::vector<int>& small_primes) {
    const auto factors = factorize(p - 1, small_primes);
    for (int g = 2; g < p; ++g) {
        bool ok = true;
        for (auto [q, exponent] : factors) {
            (void)exponent;
            long long value = 1;
            long long base = g;
            int power = (p - 1) / q;
            while (power > 0) {
                if (power & 1) {
                    value = (value * base) % p;
                }
                base = (base * base) % p;
                power >>= 1;
            }
            if (value == 1) {
                ok = false;
                break;
            }
        }
        if (ok) {
            return g;
        }
    }
    return -1;
}

bool is_prime_by_trial(long long n, const std::vector<int>& primes) {
    if (n < 2) {
        return false;
    }
    for (int q : primes) {
        if (1LL * q * q > n) {
            return true;
        }
        if (n % q == 0) {
            return n == q;
        }
    }
    return true;
}

int smallest_transition_prime(int p, int quotient_order, int target_gcd,
                              const std::vector<int>& discrete_log,
                              const std::vector<int>& primes) {
    auto works = [&](long long q) {
        if (q == p) {
            return false;
        }
        const int residue = static_cast<int>(q % p);
        if (residue == 0) {
            return false;
        }
        const int e = discrete_log[residue];
        return std::gcd(e, quotient_order) == target_gcd;
    };

    for (int q : primes) {
        if (works(q)) {
            return q;
        }
    }

    long long q = primes.back() + 1LL;
    if ((q & 1LL) == 0) {
        ++q;
    }
    for (;; q += 2) {
        if (is_prime_by_trial(q, primes) && works(q)) {
            return static_cast<int>(q);
        }
    }
}

struct PrimeSolution {
    long double log10_value = 0.0L;
    std::vector<std::pair<int, int>> factors;
};

PrimeSolution solve_prime(int p, const std::vector<int>& primes) {
    if (p == 2) {
        return {};
    }

    const int n = p - 1;
    const int g = primitive_root(p, primes);
    assert(g > 0);

    std::vector<int> discrete_log(p, -1);
    int x = 1;
    for (int e = 0; e < n; ++e) {
        discrete_log[x] = e;
        x = static_cast<int>(1LL * x * g % p);
    }

    const auto divisors = divisors_from_factors(factorize(n, primes));
    const int d_count = static_cast<int>(divisors.size());
    const long double inf = std::numeric_limits<long double>::infinity();

    std::vector<long double> dp(d_count, inf);
    std::vector<int> parent(d_count, -1);
    std::vector<std::pair<int, int>> edge(d_count, {0, 0});

    auto index_of = [&](int value) {
        return static_cast<int>(std::lower_bound(divisors.begin(), divisors.end(), value) - divisors.begin());
    };

    dp[index_of(1)] = 0.0L;
    for (int i = 0; i < d_count; ++i) {
        const int subgroup_size = divisors[i];
        if (!std::isfinite(dp[i])) {
            continue;
        }

        const int quotient_order = n / subgroup_size;
        for (int length : divisors) {
            if (length <= 1 || quotient_order % length != 0) {
                continue;
            }
            const int target_gcd = quotient_order / length;
            const int q = smallest_transition_prime(p, quotient_order, target_gcd, discrete_log, primes);
            const int next_size = subgroup_size * length;
            const int j = index_of(next_size);
            const long double candidate = dp[i] + static_cast<long double>(length - 1) * std::log10(static_cast<long double>(q));
            if (candidate + 1e-24L < dp[j]) {
                dp[j] = candidate;
                parent[j] = i;
                edge[j] = {length, q};
            }
        }
    }

    const int finish = index_of(n);
    PrimeSolution solution;
    solution.log10_value = dp[finish];

    for (int at = finish; parent[at] != -1; at = parent[at]) {
        solution.factors.push_back(edge[at]);
    }
    std::reverse(solution.factors.begin(), solution.factors.end());
    return solution;
}

std::uint64_t exact_value(const PrimeSolution& solution) {
    std::uint64_t value = 1;
    for (auto [length, q] : solution.factors) {
        for (int i = 1; i < length; ++i) {
            value *= static_cast<std::uint64_t>(q);
        }
    }
    return value;
}

std::string scientific(long double log10_value) {
    long long exponent = static_cast<long long>(std::floor(log10_value));
    long double mantissa = std::pow(10.0L, log10_value - static_cast<long double>(exponent));
    mantissa = std::round(mantissa * 100000.0L) / 100000.0L;
    if (mantissa >= 10.0L) {
        mantissa /= 10.0L;
        ++exponent;
    }

    std::ostringstream out;
    out << std::fixed << std::setprecision(5) << static_cast<double>(mantissa) << 'e' << exponent;
    return out.str();
}

struct Solver {
    std::vector<int> primes = sieve_primes(kPrimeSearchLimit);

    PrimeSolution S(int p) const {
        return solve_prime(p, primes);
    }

    long double log_T(int limit) const {
        long double total = 0.0L;
        for (int p : primes) {
            if (p >= limit) {
                break;
            }
            total += S(p).log10_value;
        }
        return total;
    }
};

void validate() {
    Solver solver;

    assert(exact_value(solver.S(2)) == 1);
    assert(exact_value(solver.S(5)) == 8);

    std::uint64_t t20 = 1;
    for (int p : solver.primes) {
        if (p >= 20) {
            break;
        }
        t20 *= exact_value(solver.S(p));
    }
    assert(t20 == 1348422598656ULL);

    assert(scientific(solver.log_T(100)) == "1.37451e123");
    std::cerr << "Validation checkpoints passed.\n";
}

}  // namespace

int main() {
    validate();

    Solver solver;
    std::cout << scientific(solver.log_T(kLimit)) << '\n';
    return 0;
}
