#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'007ULL;

u64 add_mod(u64 a, u64 b) {
    a += b;
    if (a >= kMod) {
        a -= kMod;
    }
    return a;
}

u64 mul_mod(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

u64 pow_mod(u64 base, u64 exp) {
    u64 result = 1;
    while (exp > 0) {
        if (exp & 1ULL) {
            result = mul_mod(result, base);
        }
        base = mul_mod(base, base);
        exp >>= 1ULL;
    }
    return result;
}

u64 inv_mod(u64 x) {
    return pow_mod(x, kMod - 2);
}

struct Event {
    int k;
    u64 mul;

    bool operator<(const Event& other) const {
        return k < other.k;
    }
};

std::vector<int> primes_up_to(int n) {
    std::vector<bool> composite(n + 1, false);
    std::vector<int> primes;
    primes.reserve(n / 10);

    for (int i = 2; i <= n; ++i) {
        if (!composite[i]) {
            primes.push_back(i);
            if (static_cast<i64>(i) * i <= n) {
                for (int j = i * i; j <= n; j += i) {
                    composite[j] = true;
                }
            }
        }
    }
    return primes;
}

int vp_factorial(int n, int p) {
    int e = 0;
    while (n > 0) {
        n /= p;
        e += n;
    }
    return e;
}

u64 r_from_exponents_direct(const std::vector<int>& exponents) {
    int max_e = 0;
    for (int e : exponents) {
        max_e = std::max(max_e, e);
    }

    u64 total = 0;
    for (int k = 1; k <= max_e; ++k) {
        u64 d = 1;
        for (int e : exponents) {
            d = mul_mod(d, static_cast<u64>(e / k + 1));
        }
        total = add_mod(total, d - 1);
    }
    return total;
}

std::vector<int> factor_exponents_u64(u64 n) {
    std::vector<int> exponents;
    for (u64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) {
            continue;
        }
        int c = 0;
        while (n % p == 0) {
            n /= p;
            ++c;
        }
        exponents.push_back(c);
    }
    if (n > 1) {
        exponents.push_back(1);
    }
    return exponents;
}

u64 r_of_number_direct(u64 n) {
    if (n <= 1) {
        return 0;
    }
    return r_from_exponents_direct(factor_exponents_u64(n));
}

u64 r_of_factorial_optimized(int n) {
    if (n < 2) {
        return 0;
    }

    const std::vector<int> primes = primes_up_to(n);

    std::unordered_map<int, int> count_by_exp;
    count_by_exp.reserve(primes.size() * 2);

    int max_e = 0;
    for (int p : primes) {
        const int e = vp_factorial(n, p);
        ++count_by_exp[e];
        max_e = std::max(max_e, e);
    }

    std::vector<std::pair<int, int>> groups;
    groups.reserve(count_by_exp.size());
    for (const auto& kv : count_by_exp) {
        groups.push_back(kv);
    }
    std::sort(groups.begin(), groups.end());

    u64 d1 = 1;
    for (const auto& [e, c] : groups) {
        d1 = mul_mod(d1, pow_mod(static_cast<u64>(e + 1), static_cast<u64>(c)));
    }

    std::vector<Event> events;
    events.reserve(groups.size() * 300);

    for (const auto& [e, c] : groups) {
        u64 prev_term = static_cast<u64>(e + 1);

        int l = 2;
        while (l <= e) {
            const int q = e / l;
            const int r = e / q;
            const u64 term = static_cast<u64>(q + 1);
            if (term != prev_term) {
                const u64 ratio = mul_mod(term, inv_mod(prev_term));
                events.push_back({l, pow_mod(ratio, static_cast<u64>(c))});
                prev_term = term;
            }
            l = r + 1;
        }

        if (e + 1 <= max_e && prev_term != 1) {
            const u64 ratio = inv_mod(prev_term);
            events.push_back({e + 1, pow_mod(ratio, static_cast<u64>(c))});
        }
    }

    std::sort(events.begin(), events.end());

    u64 curr_d = d1;
    u64 answer = curr_d - 1;

    std::size_t idx = 0;
    for (int k = 2; k <= max_e; ++k) {
        while (idx < events.size() && events[idx].k == k) {
            curr_d = mul_mod(curr_d, events[idx].mul);
            ++idx;
        }
        answer = add_mod(answer, curr_d - 1);
    }

    return answer;
}

std::vector<int> factorial_exponents_direct(int n) {
    const std::vector<int> primes = primes_up_to(n);
    std::vector<int> exponents;
    exponents.reserve(primes.size());
    for (int p : primes) {
        exponents.push_back(vp_factorial(n, p));
    }
    return exponents;
}

void run_validations() {
    assert(r_of_number_direct(20) == 6);
    assert(r_of_factorial_optimized(10) == 312);

    for (int n : {2, 3, 5, 8, 10, 20, 50, 100}) {
        const u64 optimized = r_of_factorial_optimized(n);
        const u64 direct = r_from_exponents_direct(factorial_exponents_direct(n));
        assert(optimized == direct);
    }
}

}  // namespace

int main() {
    run_validations();
    constexpr int kN = 10'000'000;
    std::cout << r_of_factorial_optimized(kN) << '\n';
    return 0;
}
