#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <future>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using std::uint64_t;
using u128 = __uint128_t;

namespace {

constexpr uint64_t kDefaultN = 100000000000000ULL;

uint64_t gcd_u64(uint64_t a, uint64_t b) {
    while (b) {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

uint64_t mul_mod_u64(uint64_t a, uint64_t b, uint64_t mod) {
    return static_cast<uint64_t>((static_cast<u128>(a) * b) % mod);
}

uint64_t pow_mod_u64(uint64_t a, uint64_t d, uint64_t mod) {
    uint64_t r = 1;
    while (d) {
        if (d & 1) r = mul_mod_u64(r, a, mod);
        a = mul_mod_u64(a, a, mod);
        d >>= 1;
    }
    return r;
}

bool is_prime_u64(uint64_t n) {
    if (n < 2) return false;
    static uint64_t small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (uint64_t p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    uint64_t d = n - 1, s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    auto witness = [&](uint64_t a) -> bool {
        if (a % n == 0) return false;
        uint64_t x = pow_mod_u64(a, d, n);
        if (x == 1 || x == n - 1) return false;
        for (uint64_t i = 1; i < s; ++i) {
            x = mul_mod_u64(x, x, n);
            if (x == n - 1) return false;
        }
        return true;
    };

    static uint64_t bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};
    for (uint64_t a : bases) {
        if (witness(a)) return false;
    }
    return true;
}

uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

thread_local uint64_t rng_state =
    static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count()) ^
    static_cast<uint64_t>(std::hash<std::thread::id>()(std::this_thread::get_id()));

uint64_t rand_u64(uint64_t lo, uint64_t hi) {
    uint64_t r = splitmix64(rng_state);
    return lo + (hi > lo ? (r % (hi - lo + 1)) : 0);
}

uint64_t pollard_rho(uint64_t n) {
    if ((n & 1ULL) == 0) return 2;
    if (n % 3ULL == 0) return 3;

    while (true) {
        uint64_t c = rand_u64(1, n - 1);
        uint64_t x = rand_u64(0, n - 1);
        uint64_t y = x;
        uint64_t d = 1;

        auto f = [&](uint64_t v) { return (mul_mod_u64(v, v, n) + c) % n; };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            uint64_t diff = (x > y) ? (x - y) : (y - x);
            d = gcd_u64(diff, n);
        }
        if (d != n) return d;
    }
}

void factor_u64(uint64_t n, std::vector<uint64_t>& fac) {
    if (n == 1) return;
    if (is_prime_u64(n)) {
        fac.push_back(n);
        return;
    }
    uint64_t d = pollard_rho(n);
    factor_u64(d, fac);
    factor_u64(n / d, fac);
}

struct FactorCache {
    std::unordered_map<uint64_t, std::vector<uint64_t>> mp;

    std::vector<uint64_t> get(uint64_t n) {
        auto it = mp.find(n);
        if (it != mp.end()) return it->second;

        std::vector<uint64_t> f;
        if (n > 1) {
            std::vector<uint64_t> tmp;
            factor_u64(n, tmp);
            std::sort(tmp.begin(), tmp.end());
            tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
            f = std::move(tmp);
        }
        mp.emplace(n, f);
        return f;
    }
};

int v3_u64(uint64_t x) {
    int c = 0;
    while (x % 3ULL == 0) {
        x /= 3ULL;
        ++c;
    }
    return c;
}

uint64_t sigma_prime_power(uint64_t p, int e, uint64_t p_pow) {
    u128 next_pow = static_cast<u128>(p_pow) * p;
    u128 sig = (next_pow - 1) / (p - 1);
    return static_cast<uint64_t>(sig);
}

std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct SolverA {
    uint64_t A;
    uint64_t limit;
    int maxv3;
    bool validate;

    uint64_t sum_m = 0;

    std::unordered_set<uint64_t> visited_m;
    FactorCache cache;
    std::vector<uint64_t> used;

    SolverA(uint64_t A_, uint64_t limit_, int maxv3_, bool validate_)
        : A(A_), limit(limit_), maxv3(maxv3_), validate(validate_) {
        visited_m.reserve(1 << 16);
        cache.mp.reserve(1 << 16);
        used.reserve(32);
    }

    bool is_used(uint64_t p) const {
        for (uint64_t x : used) {
            if (x == p) return true;
        }
        return false;
    }

    void dfs(uint64_t m, uint64_t num, uint64_t den, int v3sig, uint64_t sigma_m) {
        if (m > limit) return;
        if (!visited_m.insert(m).second) return;

        if (den == 1 && v3sig <= maxv3) {
            sum_m += m;

            if (validate) {
                uint64_t g = gcd_u64(m, sigma_m);
                uint64_t r = m / g;
                if (A % r != 0) {
                    std::cerr << "Validation failed: A%r!=0, A=" << A << " m=" << m << "\n";
                    std::abort();
                }
                if (v3_u64(sigma_m) != v3sig) {
                    std::cerr << "Validation failed: v3 mismatch, m=" << m << "\n";
                    std::abort();
                }
                u128 prod = static_cast<u128>(A) * sigma_m;
                if (static_cast<uint64_t>(prod % m) != 0ULL) {
                    std::cerr << "Validation failed: (A*sigma_m)%m!=0, m=" << m << "\n";
                    std::abort();
                }
            }
        }

        std::vector<uint64_t> cand = cache.get(num);
        cand.push_back(2);
        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

        for (uint64_t p : cand) {
            if (p == 3) continue;
            if (is_used(p)) continue;

            uint64_t p_pow = p;
            for (int e = 1;; ++e) {
                if (static_cast<u128>(m) * p_pow > limit) break;

                uint64_t sig = sigma_prime_power(p, e, p_pow);
                int v3f = v3_u64(sig);
                if (v3sig + v3f <= maxv3) {
                    uint64_t denom_factor = p_pow;
                    uint64_t sig_factor = sig;

                    uint64_t num1 = num;
                    uint64_t den1 = den;

                    uint64_t g1 = gcd_u64(num1, denom_factor);
                    num1 /= g1;
                    denom_factor /= g1;

                    uint64_t g2 = gcd_u64(sig_factor, den1);
                    sig_factor /= g2;
                    den1 /= g2;

                    u128 new_num = static_cast<u128>(num1) * sig_factor;
                    u128 new_den = static_cast<u128>(den1) * denom_factor;

                    uint64_t new_m = static_cast<uint64_t>(static_cast<u128>(m) * p_pow);
                    uint64_t new_sigma_m = static_cast<uint64_t>(static_cast<u128>(sigma_m) * sig);

                    used.push_back(p);
                    dfs(new_m, static_cast<uint64_t>(new_num), static_cast<uint64_t>(new_den),
                        v3sig + v3f, new_sigma_m);
                    used.pop_back();
                }

                if (p_pow > limit / p) break;
                p_pow *= p;
            }
        }
    }

    void run() {
        dfs(1, A, 1, 0, 1);
    }
};

u128 computeT(uint64_t N, bool validate, bool multithread) {
    std::vector<std::pair<int, uint64_t>> tasks;
    u128 pow3 = 1;
    for (int a = 1;; ++a) {
        pow3 *= 3;
        if (pow3 > N) break;
        tasks.push_back({a, static_cast<uint64_t>(pow3)});
    }

    auto solve_one = [&](int a, uint64_t p3) -> u128 {
        uint64_t limit = N / p3;
        u128 t = static_cast<u128>(p3) * 3 - 1;
        uint64_t A = static_cast<uint64_t>(t / 2);
        SolverA solver(A, limit, a - 1, validate);
        solver.run();
        return static_cast<u128>(p3) * solver.sum_m;
    };

    u128 total = 0;
    if (!multithread || tasks.size() <= 1) {
        for (const auto& task : tasks) total += solve_one(task.first, task.second);
        return total;
    }

    std::vector<std::future<u128>> fut;
    fut.reserve(tasks.size());
    for (const auto& task : tasks) {
        fut.push_back(std::async(std::launch::async, solve_one, task.first, task.second));
    }
    for (auto& f : fut) total += f.get();
    return total;
}

}  // namespace

int main(int argc, char** argv) {
    uint64_t N = kDefaultN;
    bool validate = true;
    bool multithread = true;

    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if (s.rfind("--N=", 0) == 0) {
            N = std::stoull(s.substr(4));
        } else if (s.rfind("--validate=", 0) == 0) {
            validate = (std::stoi(s.substr(11)) != 0);
        } else if (s.rfind("--mt=", 0) == 0) {
            multithread = (std::stoi(s.substr(5)) != 0);
        }
    }

    if (validate) {
        u128 t100 = computeT(100ULL, false, false);
        if (static_cast<uint64_t>(t100) != 270ULL) {
            std::cerr << "Checkpoint failed: T(100) expected 270, got "
                      << static_cast<uint64_t>(t100) << "\n";
            return 1;
        }
        u128 t1e6 = computeT(1000000ULL, false, false);
        if (static_cast<uint64_t>(t1e6) != 26089287ULL) {
            std::cerr << "Checkpoint failed: T(1e6) expected 26089287, got "
                      << static_cast<uint64_t>(t1e6) << "\n";
            return 1;
        }
    }

    u128 ans = computeT(N, validate, multithread);
    std::cout << to_string_u128(ans) << "\n";
    return 0;
}
