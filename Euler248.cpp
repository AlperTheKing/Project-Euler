#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kTargetPhi = 6'227'020'800ULL;  // 13!
constexpr u64 kKnownFirstSolution = 6'227'180'929ULL;
constexpr std::size_t kDefaultIndex = 150'000;

struct SeedTask {
    u64 rem_phi = 1;
    std::size_t next_index = 0;
    u128 current_n = 1;
};

u64 pow_u64(u64 base, int exp) {
    u64 out = 1;
    for (int i = 0; i < exp; ++i) out *= base;
    return out;
}

u64 mul_mod_u64(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod_u64(u64 base, u64 exp, u64 mod) {
    u64 out = 1 % mod;
    u64 cur = base % mod;
    while (exp > 0) {
        if (exp & 1ULL) out = mul_mod_u64(out, cur, mod);
        cur = mul_mod_u64(cur, cur, mod);
        exp >>= 1ULL;
    }
    return out;
}

bool is_prime_u64(u64 n) {
    if (n < 2) return false;
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0) {
        d >>= 1ULL;
        ++s;
    }

    // Deterministic bases for 64-bit integers.
    constexpr u64 kBases[] = {
        2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL,
    };

    for (u64 a : kBases) {
        if (a % n == 0) continue;
        u64 x = pow_mod_u64(a, d, n);
        if (x == 1 || x == n - 1) continue;

        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod_u64(x, x, n);
            if (x == n - 1) {
                witness = false;
                break;
            }
        }
        if (witness) return false;
    }

    return true;
}

std::vector<std::pair<u64, int>> factorize_u64(u64 n) {
    std::vector<std::pair<u64, int>> fac;
    for (u64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) continue;
        int e = 0;
        while (n % p == 0) {
            n /= p;
            ++e;
        }
        fac.push_back({p, e});
    }
    if (n > 1) fac.push_back({n, 1});
    return fac;
}

void build_divisors_dfs(const std::vector<std::pair<u64, int>>& fac,
                        std::size_t idx,
                        u64 cur,
                        std::vector<u64>& out) {
    if (idx == fac.size()) {
        out.push_back(cur);
        return;
    }

    const u64 p = fac[idx].first;
    const int e = fac[idx].second;
    u64 value = 1;
    for (int i = 0; i <= e; ++i) {
        build_divisors_dfs(fac, idx + 1, cur * value, out);
        value *= p;
    }
}

std::vector<u64> divisors_of(u64 n) {
    const auto fac = factorize_u64(n);
    std::vector<u64> out;
    out.reserve(2048);
    build_divisors_dfs(fac, 0, 1, out);
    std::sort(out.begin(), out.end());
    return out;
}

std::vector<u64> candidate_primes_for_inverse_totient(u64 phi_target) {
    std::vector<u64> candidates;
    const std::vector<u64> divs = divisors_of(phi_target);
    candidates.reserve(divs.size());

    for (u64 d : divs) {
        const u64 p = d + 1;
        if (is_prime_u64(p)) candidates.push_back(p);
    }

    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    return candidates;
}

void dfs_inverse_totient(u64 rem_phi,
                         std::size_t start_index,
                         u128 current_n,
                         const std::vector<u64>& candidates,
                         std::vector<u128>& out) {
    if (rem_phi == 1) {
        out.push_back(current_n);
        return;
    }

    for (std::size_t i = start_index; i < candidates.size(); ++i) {
        const u64 p = candidates[i];
        u64 phi_piece = p - 1;
        if (rem_phi % phi_piece != 0) continue;

        u128 n_piece = static_cast<u128>(p);
        while (rem_phi % phi_piece == 0) {
            dfs_inverse_totient(rem_phi / phi_piece, i + 1, current_n * n_piece, candidates, out);

            if (phi_piece > rem_phi / p) break;
            phi_piece *= p;
            n_piece *= static_cast<u128>(p);
        }
    }
}

std::vector<SeedTask> build_seed_tasks(u64 phi_target, const std::vector<u64>& candidates) {
    std::vector<SeedTask> tasks;
    tasks.reserve(candidates.size() * 2);

    for (std::size_t i = 0; i < candidates.size(); ++i) {
        const u64 p = candidates[i];
        u64 phi_piece = p - 1;
        u128 n_piece = static_cast<u128>(p);

        while (phi_target % phi_piece == 0) {
            tasks.push_back(SeedTask{phi_target / phi_piece, i + 1, n_piece});
            if (phi_piece > phi_target / p) break;
            phi_piece *= p;
            n_piece *= static_cast<u128>(p);
        }
    }

    return tasks;
}

std::vector<u128> enumerate_inverse_totient(u64 phi_target,
                                            const std::vector<u64>& candidates,
                                            int threads) {
    if (threads < 1) threads = 1;

    std::vector<u128> out;
    if (threads == 1) {
        dfs_inverse_totient(phi_target, 0, 1, candidates, out);
        if (phi_target == 1) out.push_back(1);
        return out;
    }

    const std::vector<SeedTask> tasks = build_seed_tasks(phi_target, candidates);
    std::atomic<std::size_t> next_task{0};
    std::vector<std::vector<u128>> partial(static_cast<std::size_t>(threads));
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            std::vector<u128>& local = partial[static_cast<std::size_t>(t)];
            while (true) {
                const std::size_t id = next_task.fetch_add(1, std::memory_order_relaxed);
                if (id >= tasks.size()) break;
                const SeedTask& task = tasks[id];
                dfs_inverse_totient(task.rem_phi, task.next_index, task.current_n, candidates, local);
            }
        });
    }

    for (std::thread& th : pool) th.join();

    std::size_t total = (phi_target == 1) ? 1 : 0;
    for (const auto& vec : partial) total += vec.size();
    out.reserve(total);

    if (phi_target == 1) out.push_back(1);
    for (auto& vec : partial) {
        out.insert(out.end(), vec.begin(), vec.end());
    }
    return out;
}

void normalize_solutions(std::vector<u128>& solutions) {
    std::sort(solutions.begin(), solutions.end());
    solutions.erase(std::unique(solutions.begin(), solutions.end()), solutions.end());
}

u64 phi_u64(u64 n) {
    if (n == 0) return 0;
    u64 result = n;

    for (u64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) continue;
        while (n % p == 0) n /= p;
        result = result / p * (p - 1);
    }
    if (n > 1) result = result / n * (n - 1);
    return result;
}

u64 phi_from_candidate_factorization(u128 n, const std::vector<u64>& candidates) {
    u128 result = 1;
    u128 rem = n;

    for (u64 p : candidates) {
        if (rem % p != 0) continue;

        int exp = 0;
        while (rem % p == 0) {
            rem /= p;
            ++exp;
        }

        u128 term = static_cast<u128>(p - 1);
        for (int i = 1; i < exp; ++i) term *= static_cast<u128>(p);
        result *= term;
    }

    if (rem != 1) return 0;
    if (result > std::numeric_limits<u64>::max()) return 0;
    return static_cast<u64>(result);
}

std::string to_string_u128(u128 value) {
    if (value == 0) return "0";
    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

bool validate_small_case() {
    constexpr u64 small_target = 120;
    constexpr u64 brute_limit = 200'000;

    const std::vector<u64> candidates = candidate_primes_for_inverse_totient(small_target);
    std::vector<u128> by_solver = enumerate_inverse_totient(small_target, candidates, 1);
    normalize_solutions(by_solver);

    std::vector<u128> by_bruteforce;
    for (u64 n = 1; n <= brute_limit; ++n) {
        if (phi_u64(n) == small_target) by_bruteforce.push_back(static_cast<u128>(n));
    }

    if (by_solver != by_bruteforce) {
        std::cerr << "Validation failed for phi=120 (solver vs brute force mismatch).\n";
        std::cerr << "solver size=" << by_solver.size() << ", brute size=" << by_bruteforce.size() << "\n";
        return false;
    }

    return true;
}

bool validate_main_problem(const std::vector<u64>& candidates,
                           const std::vector<u128>& threaded_solutions,
                           int threads) {
    if (threaded_solutions.empty()) {
        std::cerr << "Validation failed: no solutions found for phi=13!.\n";
        return false;
    }

    if (threaded_solutions.front() != static_cast<u128>(kKnownFirstSolution)) {
        std::cerr << "Validation failed: first solution mismatch. got="
                  << to_string_u128(threaded_solutions.front())
                  << " expected=" << kKnownFirstSolution << "\n";
        return false;
    }

    for (u128 n : threaded_solutions) {
        const u64 phi_value = phi_from_candidate_factorization(n, candidates);
        if (phi_value != kTargetPhi) {
            std::cerr << "Validation failed: phi(n) mismatch for n=" << to_string_u128(n)
                      << ", got phi=" << phi_value << "\n";
            return false;
        }
    }

    if (threads > 1) {
        std::vector<u128> single = enumerate_inverse_totient(kTargetPhi, candidates, 1);
        normalize_solutions(single);
        if (single != threaded_solutions) {
            std::cerr << "Validation failed: single-thread and multi-thread outputs differ.\n";
            std::cerr << "single size=" << single.size() << ", multi size=" << threaded_solutions.size() << "\n";
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::size_t index = kDefaultIndex;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 4;

    if (argc > 1) {
        const long long parsed = std::atoll(argv[1]);
        if (parsed > 0) index = static_cast<std::size_t>(parsed);
    }
    if (argc > 2) {
        const int parsed_threads = std::atoi(argv[2]);
        if (parsed_threads > 0) threads = parsed_threads;
    }

    if (!validate_small_case()) return 1;

    const std::vector<u64> candidates = candidate_primes_for_inverse_totient(kTargetPhi);
    std::vector<u128> solutions = enumerate_inverse_totient(kTargetPhi, candidates, threads);
    normalize_solutions(solutions);

    if (!validate_main_problem(candidates, solutions, threads)) return 1;

    if (index == 0 || index > solutions.size()) {
        std::cerr << "Requested index out of range: " << index
                  << " (number of solutions=" << solutions.size() << ")\n";
        return 1;
    }

    const u128 answer = solutions[index - 1];
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
