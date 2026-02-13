#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        s.push_back(static_cast<char>('0' + static_cast<int>(x % 10)));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::vector<int> sieve_primes(int limit) {
    std::vector<bool> composite(limit + 1, false);
    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (!composite[i]) {
            primes.push_back(i);
            if (static_cast<u64>(i) * static_cast<u64>(i) <= static_cast<u64>(limit)) {
                for (int j = i * i; j <= limit; j += i) {
                    composite[j] = true;
                }
            }
        }
    }
    return primes;
}

u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) {
            r = static_cast<u64>((static_cast<u128>(r) * a) % mod);
        }
        a = static_cast<u64>((static_cast<u128>(a) * a) % mod);
        e >>= 1ULL;
    }
    return r;
}

u64 mod_inv_prime(u64 a, u64 p) {
    return mod_pow(a, p - 2, p);
}

u128 solve(u64 N) {
    const std::vector<int> primes = sieve_primes(10'000);

    std::vector<u64> residues = {0};
    u128 modulus = 1;

    u128 pass_count = N;
    u128 answer = 0;

    for (int p_int : primes) {
        const u64 p = static_cast<u64>(p_int);
        const u128 next_modulus = modulus * p;

        const u64 inv = mod_inv_prime(static_cast<u64>(modulus % p), p);

        std::vector<u64> next_residues;
        next_residues.reserve(residues.size() * (p / 7 + 1));

        for (u64 r : residues) {
            const u64 rp = r % p;
            for (u64 a = 0; a < p; a += 7) {
                const u64 t = static_cast<u64>((static_cast<u128>((a + p - rp) % p) * inv) % p);

                if (modulus > N && t != 0) {
                    continue;
                }

                const u128 x = static_cast<u128>(r) + modulus * t;
                if (x <= N) {
                    next_residues.push_back(static_cast<u64>(x));
                }
            }
        }

        u128 next_count = 0;
        if (next_modulus <= N) {
            for (u64 r : next_residues) {
                if (r == 0) {
                    next_count += N / next_modulus;
                } else {
                    next_count += (N - r) / next_modulus + 1;
                }
            }
        } else {
            for (u64 r : next_residues) {
                if (r != 0) {
                    ++next_count;
                }
            }
        }

        assert(next_count <= pass_count);
        const u128 fail_here = pass_count - next_count;
        answer += fail_here * p;

        if (next_count == 0) {
            return answer;
        }

        residues.swap(next_residues);
        modulus = next_modulus;
        pass_count = next_count;
    }

    assert(false);
    return 0;
}

u128 brute(u64 N) {
    const std::vector<int> primes = sieve_primes(200'000);

    u128 total = 0;
    for (u64 n = 1; n <= N; ++n) {
        bool found = false;
        for (int p : primes) {
            if ((n % static_cast<u64>(p)) % 7ULL != 0ULL) {
                total += static_cast<u64>(p);
                found = true;
                break;
            }
        }
        assert(found);
    }
    return total;
}

void run_validations() {
    assert(solve(1470) == 4293);

    for (u64 n : {10ULL, 100ULL, 1000ULL, 10'000ULL}) {
        assert(solve(n) == brute(n));
    }
}

}  // namespace

int main() {
    run_validations();
    constexpr u64 kN = 100'000'000'000'000'000ULL;
    std::cout << to_string_u128(solve(kN)) << '\n';
    return 0;
}
