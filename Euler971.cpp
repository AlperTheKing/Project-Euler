#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1ULL) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp >>= 1ULL;
    }
    return result;
}

std::vector<int> primes_up_to(int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }

    for (int i = 2; static_cast<long long>(i) * i <= n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = false;
        }
    }

    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

int count_cycle_nodes(const std::array<int, 5>& next) {
    std::array<int, 5> state{};  // 0=unseen, 1=done, 2=in-stack
    int cycle_nodes = 0;

    for (int start = 0; start < 5; ++start) {
        if (state[static_cast<std::size_t>(start)] != 0) {
            continue;
        }

        std::array<int, 5> stack{};
        std::array<int, 5> pos{};
        pos.fill(-1);
        int top = 0;

        int cur = start;
        while (state[static_cast<std::size_t>(cur)] == 0) {
            state[static_cast<std::size_t>(cur)] = 2;
            pos[static_cast<std::size_t>(cur)] = top;
            stack[static_cast<std::size_t>(top++)] = cur;
            cur = next[static_cast<std::size_t>(cur)];
        }

        if (state[static_cast<std::size_t>(cur)] == 2) {
            cycle_nodes += top - pos[static_cast<std::size_t>(cur)];
        }

        for (int i = 0; i < top; ++i) {
            state[static_cast<std::size_t>(stack[static_cast<std::size_t>(i)])] = 1;
        }
    }

    return cycle_nodes;
}

u64 count_periodic_points_for_prime(int p) {
    const u64 t = static_cast<u64>(p - 1) / 5ULL;

    u64 omega = 1;
    for (u64 a = 2; a < static_cast<u64>(p); ++a) {
        omega = mod_pow(a, t, static_cast<u64>(p));
        if (omega != 1) {
            break;
        }
    }

    std::array<u64, 5> pow_omega{};
    pow_omega[0] = 1;
    for (int i = 1; i < 5; ++i) {
        pow_omega[static_cast<std::size_t>(i)] =
            (pow_omega[static_cast<std::size_t>(i - 1)] * omega) % static_cast<u64>(p);
    }

    std::array<int, 5> next{};
    for (int r = 0; r < 5; ++r) {
        const u64 value = mod_pow((1 + pow_omega[static_cast<std::size_t>(r)]) % static_cast<u64>(p),
                                  t, static_cast<u64>(p));
        int shift = -1;
        for (int i = 0; i < 5; ++i) {
            if (pow_omega[static_cast<std::size_t>(i)] == value) {
                shift = i;
                break;
            }
        }
        if (shift < 0) {
            return 0ULL;
        }
        next[static_cast<std::size_t>(r)] = (r + shift) % 5;
    }

    const int cycle_nodes = count_cycle_nodes(next);
    return 1ULL + t * static_cast<u64>(cycle_nodes);
}

u64 solve(int limit) {
    const std::vector<int> primes = primes_up_to(limit);

    u64 total = 0;
    for (int p : primes) {
        if (p % 5 != 1) {
            continue;
        }
        total += count_periodic_points_for_prime(p);
    }
    return total;
}

u64 brute_count_periodic_points_for_prime(int p) {
    const int k = (p + 4) / 5;
    std::vector<int> next(static_cast<std::size_t>(p), 0);
    for (int x = 0; x < p; ++x) {
        next[static_cast<std::size_t>(x)] =
            static_cast<int>((mod_pow(static_cast<u64>(x), static_cast<u64>(k),
                                      static_cast<u64>(p)) +
                              static_cast<u64>(x)) %
                             static_cast<u64>(p));
    }

    std::vector<int> state(static_cast<std::size_t>(p), 0);  // 0=unseen,1=done,2=in-stack
    int periodic = 0;

    for (int start = 0; start < p; ++start) {
        if (state[static_cast<std::size_t>(start)] != 0) {
            continue;
        }

        std::vector<int> stack;
        stack.reserve(64);
        std::vector<int> pos(static_cast<std::size_t>(p), -1);

        int cur = start;
        while (state[static_cast<std::size_t>(cur)] == 0) {
            state[static_cast<std::size_t>(cur)] = 2;
            pos[static_cast<std::size_t>(cur)] = static_cast<int>(stack.size());
            stack.push_back(cur);
            cur = next[static_cast<std::size_t>(cur)];
        }

        if (state[static_cast<std::size_t>(cur)] == 2) {
            periodic += static_cast<int>(stack.size()) - pos[static_cast<std::size_t>(cur)];
        }

        for (int v : stack) {
            state[static_cast<std::size_t>(v)] = 1;
        }
    }

    return static_cast<u64>(periodic);
}

bool run_checkpoints() {
    if (count_periodic_points_for_prime(11) != 7ULL) {
        std::cerr << "Checkpoint failed: C(11) != 7" << '\n';
        return false;
    }

    for (int p : {11, 31, 41, 61}) {
        const u64 fast = count_periodic_points_for_prime(p);
        const u64 brute = brute_count_periodic_points_for_prime(p);
        if (fast != brute) {
            std::cerr << "Checkpoint failed for p=" << p << ": fast=" << fast
                      << ", brute=" << brute << '\n';
            return false;
        }
    }

    if (solve(100) != 127ULL) {
        std::cerr << "Checkpoint failed: S(100) != 127" << '\n';
        return false;
    }

    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    std::cout << solve(100'000'000) << '\n';
    return 0;
}
