#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;
using boost::multiprecision::cpp_int;

constexpr u64 kBase = 10'000'000'000ULL;

cpp_int mod_pow_cpp_int(cpp_int base, u64 exp, const cpp_int& mod) {
    base %= mod;
    cpp_int result = 1;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp >>= 1ULL;
    }
    return result;
}

u64 A_block(const u64 n) {
    const u64 n_plus_9 = n + 9ULL;

    std::vector<u64> p3(1, 1ULL);
    while (p3.back() <= n_plus_9 / 3ULL) {
        p3.push_back(p3.back() * 3ULL);
    }

    int limit = 0;
    for (int i = 1; i < static_cast<int>(p3.size()); ++i) {
        if (p3[static_cast<std::size_t>(i)] <= n_plus_9) {
            limit = i;
        }
    }

    if (limit == 0) {
        return 0ULL;
    }

    u64 sum_q_mod = 0ULL;
    std::vector<u64> remainders(static_cast<std::size_t>(limit + 1), 0ULL);

    for (int i = 1; i <= limit; ++i) {
        const u64 d = p3[static_cast<std::size_t>(i)];
        const u64 m = n_plus_9 - d;

        const cpp_int modulus = cpp_int(d) * cpp_int(kBase);
        const cpp_int s = mod_pow_cpp_int(cpp_int(10), m, modulus);

        const cpp_int q = s / d;
        const cpp_int r = s % d;
        const u64 q_mod = (q % kBase).convert_to<u64>();
        const u64 r_u64 = r.convert_to<u64>();

        sum_q_mod += q_mod;
        if (sum_q_mod >= kBase) {
            sum_q_mod %= kBase;
        }

        remainders[static_cast<std::size_t>(i)] = r_u64;
    }

    const u64 denom = p3[static_cast<std::size_t>(limit)];
    u128 numer = 0;
    for (int i = 1; i <= limit; ++i) {
        const u64 scale = p3[static_cast<std::size_t>(limit - i)];
        numer += static_cast<u128>(remainders[static_cast<std::size_t>(i)]) * static_cast<u128>(scale);
    }

    const u64 carry = static_cast<u64>(numer / static_cast<u128>(denom));
    return (sum_q_mod + carry) % kBase;
}

}  // namespace

int main() {
    assert(A_block(100ULL) == 4'938'271'604ULL);
    assert(A_block(100'000'000ULL) == 2'584'642'393ULL);

    const u64 answer = A_block(10'000'000'000'000'000ULL);
    std::cout << std::setw(10) << std::setfill('0') << answer << '\n';
    return 0;
}
