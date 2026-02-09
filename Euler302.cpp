#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 limit = 1000000000000000000ULL;
    unsigned threads = std::thread::hardware_concurrency();
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<u64>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg, const std::string& prefix, unsigned& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    unsigned parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10U + static_cast<unsigned>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.limit < 2) {
        return false;
    }
    if (options.threads == 0) {
        options.threads = 1;
    }
    return true;
}

int cube_root_floor(const u64 n) {
    long double x = std::cbrt(static_cast<long double>(n));
    u64 r = static_cast<u64>(x);
    while ((static_cast<u128>(r + 1) * (r + 1) * (r + 1)) <= n) {
        ++r;
    }
    while ((static_cast<u128>(r) * r * r) > n) {
        --r;
    }
    return static_cast<int>(r);
}

struct PrimeData {
    std::vector<int> primes;
    std::vector<int> spf;
    std::vector<std::vector<std::pair<int, int>>> pminus1_fact;
};

PrimeData build_prime_data(const int pmax) {
    PrimeData data;
    data.spf.assign(static_cast<std::size_t>(pmax + 1), 0);

    for (int i = 2; i <= pmax; ++i) {
        if (data.spf[static_cast<std::size_t>(i)] != 0) {
            continue;
        }
        data.spf[static_cast<std::size_t>(i)] = i;
        data.primes.push_back(i);
        if (static_cast<u64>(i) * i > static_cast<u64>(pmax)) {
            continue;
        }
        for (int j = i * i; j <= pmax; j += i) {
            if (data.spf[static_cast<std::size_t>(j)] == 0) {
                data.spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }

    data.pminus1_fact.resize(data.primes.size());

    for (std::size_t idx = 0; idx < data.primes.size(); ++idx) {
        int x = data.primes[idx] - 1;
        auto& out = data.pminus1_fact[idx];
        while (x > 1) {
            const int p = data.spf[static_cast<std::size_t>(x)];
            int e = 0;
            while (x % p == 0) {
                x /= p;
                ++e;
            }
            out.push_back({p, e});
        }
    }

    return data;
}

struct Worker {
    const PrimeData& data;
    u64 limit;

    std::vector<int> phi_exp;
    std::vector<int> active;
    u64 answer = 0;

    explicit Worker(const PrimeData& d, const u64 lim)
        : data(d), limit(lim), phi_exp(static_cast<std::size_t>(d.spf.size()), 0) {
        active.reserve(128);
    }

    inline void add_exp(const int p, const int delta, std::vector<std::pair<int, int>>& mods) {
        if (phi_exp[static_cast<std::size_t>(p)] == 0) {
            active.push_back(p);
        }
        phi_exp[static_cast<std::size_t>(p)] += delta;
        mods.push_back({p, delta});
    }

    inline void undo_mods(const std::vector<std::pair<int, int>>& mods) {
        for (auto it = mods.rbegin(); it != mods.rend(); ++it) {
            phi_exp[static_cast<std::size_t>(it->first)] -= it->second;
        }
        while (!active.empty() && phi_exp[static_cast<std::size_t>(active.back())] == 0) {
            active.pop_back();
        }
    }

    bool phi_is_achilles() const {
        int g = 0;
        bool any = false;
        for (int p : active) {
            const int e = phi_exp[static_cast<std::size_t>(p)];
            if (e == 0) {
                continue;
            }
            any = true;
            if (e < 2) {
                return false;
            }
            g = (g == 0) ? e : std::gcd(g, e);
        }
        return any && g == 1;
    }

    void dfs(const int start_idx, const u128 n_cur, const int gcd_exp, const int factor_count) {
        if (factor_count >= 2 && gcd_exp == 1 && phi_is_achilles()) {
            ++answer;
        }

        for (int i = start_idx; i < static_cast<int>(data.primes.size()); ++i) {
            const int p = data.primes[static_cast<std::size_t>(i)];
            if (static_cast<u128>(p) * p > static_cast<u128>(limit) / n_cur) {
                break;
            }

            u128 n_next = n_cur;
            int e = 0;
            while (n_next <= static_cast<u128>(limit) / p) {
                n_next *= p;
                ++e;
                if (e < 2) {
                    continue;
                }

                std::vector<std::pair<int, int>> mods;
                mods.reserve(8);

                add_exp(p, e - 1, mods);
                for (const auto& [q, qe] : data.pminus1_fact[static_cast<std::size_t>(i)]) {
                    add_exp(q, qe, mods);
                }

                const int gcd_next = (factor_count == 0) ? e : std::gcd(gcd_exp, e);
                dfs(i + 1, n_next, gcd_next, factor_count + 1);

                undo_mods(mods);
            }
        }
    }

    void run_root_index(const int i) {
        const int p = data.primes[static_cast<std::size_t>(i)];

        u128 n_next = 1;
        int e = 0;
        while (n_next <= static_cast<u128>(limit) / p) {
            n_next *= p;
            ++e;
            if (e < 2) {
                continue;
            }

            std::vector<std::pair<int, int>> mods;
            mods.reserve(8);
            add_exp(p, e - 1, mods);
            for (const auto& [q, qe] : data.pminus1_fact[static_cast<std::size_t>(i)]) {
                add_exp(q, qe, mods);
            }

            dfs(i + 1, n_next, e, 1);
            undo_mods(mods);
        }
    }
};

u64 count_strong_achilles(const u64 limit, unsigned threads) {
    const int pmax = cube_root_floor(limit);
    const PrimeData data = build_prime_data(pmax);

    if (threads <= 1 || data.primes.size() < 1000) {
        Worker worker(data, limit);
        for (int i = 0; i < static_cast<int>(data.primes.size()); ++i) {
            const int p = data.primes[static_cast<std::size_t>(i)];
            if (static_cast<u128>(p) * p > limit) {
                break;
            }
            worker.run_root_index(i);
        }
        return worker.answer;
    }

    const unsigned use_threads = std::min<unsigned>(threads, static_cast<unsigned>(data.primes.size()));
    std::vector<std::thread> pool;
    std::vector<u64> partial(use_threads, 0);
    pool.reserve(use_threads);

    for (unsigned t = 0; t < use_threads; ++t) {
        pool.emplace_back([&, t]() {
            Worker worker(data, limit);
            for (int i = static_cast<int>(t); i < static_cast<int>(data.primes.size()); i += static_cast<int>(use_threads)) {
                const int p = data.primes[static_cast<std::size_t>(i)];
                if (static_cast<u128>(p) * p > limit) {
                    break;
                }
                worker.run_root_index(i);
            }
            partial[static_cast<std::size_t>(t)] = worker.answer;
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    u64 total = 0;
    for (u64 v : partial) {
        total += v;
    }
    return total;
}

bool run_checkpoints() {
    if (count_strong_achilles(10000ULL, 1) != 7ULL) {
        std::cerr << "Checkpoint failed for 10^4" << '\n';
        return false;
    }
    if (count_strong_achilles(100000000ULL, 1) != 656ULL) {
        std::cerr << "Checkpoint failed for 10^8" << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const u64 answer = count_strong_achilles(options.limit, options.threads);
    std::cout << answer << '\n';
    return 0;
}
