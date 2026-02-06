#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'033ULL;
constexpr int kDefaultN = 100;

u64 add_mod(u64 a, u64 b) {
    a += b;
    if (a >= kMod) {
        a -= kMod;
    }
    return a;
}

u64 sub_mod(u64 a, u64 b) {
    return (a >= b) ? (a - b) : (a + kMod - b);
}

u64 mul_mod(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % kMod);
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 300'000) {
        return 1;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    threads = std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
    return threads;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    unsigned parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<unsigned>(c - '0');
    }

    value = parsed;
    return true;
}

bool parse_int_after_prefix(const std::string& arg, const char* prefix, int& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    int parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + (c - '0');
    }

    value = parsed;
    return true;
}

u64 brute_connected_bipartite(int m, int n) {
    const int edges = m * n;
    if (edges < 0 || edges > 20) {
        return 0;
    }

    const int vertices = m + n;
    const u64 masks = 1ULL << edges;
    std::vector<int> queue(vertices);
    std::vector<char> seen(vertices, 0);

    u64 connected_count = 0;

    for (u64 mask = 0; mask < masks; ++mask) {
        std::fill(seen.begin(), seen.end(), 0);

        int head = 0;
        int tail = 0;
        seen[0] = 1;
        queue[tail++] = 0;

        while (head < tail) {
            const int u = queue[head++];
            if (u < m) {
                const int row = u;
                for (int col = 0; col < n; ++col) {
                    const int edge_index = row * n + col;
                    if (((mask >> edge_index) & 1ULL) == 0ULL) {
                        continue;
                    }
                    const int v = m + col;
                    if (!seen[v]) {
                        seen[v] = 1;
                        queue[tail++] = v;
                    }
                }
            } else {
                const int col = u - m;
                for (int row = 0; row < m; ++row) {
                    const int edge_index = row * n + col;
                    if (((mask >> edge_index) & 1ULL) == 0ULL) {
                        continue;
                    }
                    const int v = row;
                    if (!seen[v]) {
                        seen[v] = 1;
                        queue[tail++] = v;
                    }
                }
            }
        }

        bool connected = true;
        for (const char s : seen) {
            if (!s) {
                connected = false;
                break;
            }
        }

        if (connected) {
            ++connected_count;
        }
    }

    return connected_count;
}

struct Options {
    int n = kDefaultN;
    bool allow_multithreading = true;
    unsigned requested_threads = 0;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        unsigned threads = 0;
        if (parse_unsigned_after_prefix(arg, "--threads=", threads)) {
            options.requested_threads = threads;
            continue;
        }

        int n = 0;
        if (parse_int_after_prefix(arg, "--n=", n)) {
            if (n < 1) {
                std::cerr << "Invalid --n value: " << n << '\n';
                return false;
            }
            options.n = n;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

class Euler434Solver {
public:
    Euler434Solver(int n, bool allow_multithreading, unsigned requested_threads)
        : n_(n),
          allow_multithreading_(allow_multithreading),
          requested_threads_(requested_threads),
          binom_(static_cast<std::size_t>(n_) + 1,
                 std::vector<u64>(static_cast<std::size_t>(n_) + 1, 0)),
          pow2_(static_cast<std::size_t>(n_ * n_) + 1, 0),
          r_(static_cast<std::size_t>(n_) + 1,
             std::vector<u64>(static_cast<std::size_t>(n_) + 1, 0)) {
        build_binomials();
        build_powers_of_two();
    }

    void solve() {
        if (n_ >= 1) {
            r_[1][0] = 1; // Connected component with one fixed row vertex and no column.
        }

        for (int m = 1; m <= n_; ++m) {
            const std::vector<u64> base_prev = compute_base_from_previous_rows(m);

            for (int n = 1; n <= n_; ++n) {
                u64 same_row = 0;
                for (int j = 0; j < n; ++j) {
                    same_row = add_mod(same_row,
                                       mul_mod(binom_[static_cast<std::size_t>(n)]
                                                     [static_cast<std::size_t>(j)],
                                               r_[static_cast<std::size_t>(m)]
                                                 [static_cast<std::size_t>(j)]));
                }

                u64 value = pow2_[static_cast<std::size_t>(m * n)];
                value = sub_mod(value, base_prev[static_cast<std::size_t>(n)]);
                value = sub_mod(value, same_row);
                r_[static_cast<std::size_t>(m)][static_cast<std::size_t>(n)] = value;
            }
        }
    }

    u64 r_value(int m, int n) const {
        if (m < 0 || m > n_ || n < 0 || n > n_) {
            return 0;
        }
        return r_[static_cast<std::size_t>(m)][static_cast<std::size_t>(n)];
    }

    u64 s_value(int n) const {
        if (n < 1) {
            return 0;
        }
        n = std::min(n, n_);

        u64 total = 0;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= n; ++j) {
                total = add_mod(total,
                                r_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
            }
        }
        return total;
    }

    bool run_checkpoints() const {
        struct Checkpoint {
            int m;
            int n;
            u64 expected;
            const char* label;
        };

        constexpr Checkpoint checks[] = {
            {2, 3, 19, "R(2,3)"},
            {5, 5, 23'679'901ULL, "R(5,5)"},
        };

        for (const Checkpoint& checkpoint : checks) {
            if (checkpoint.m > n_ || checkpoint.n > n_) {
                continue;
            }

            const u64 got = r_value(checkpoint.m, checkpoint.n);
            if (got != checkpoint.expected) {
                std::cerr << "Checkpoint failed for " << checkpoint.label
                          << ": expected " << checkpoint.expected
                          << ", got " << got << '\n';
                return false;
            }

            std::cout << checkpoint.label << " = " << got << " (ok)" << '\n';
        }

        if (n_ >= 5) {
            constexpr u64 expected_s5 = 25'021'721ULL;
            const u64 got_s5 = s_value(5);
            if (got_s5 != expected_s5) {
                std::cerr << "Checkpoint failed for S(5): expected " << expected_s5
                          << ", got " << got_s5 << '\n';
                return false;
            }
            std::cout << "S(5) = " << got_s5 << " (ok)" << '\n';
        }

        struct BruteCheckpoint {
            int m;
            int n;
        };

        constexpr BruteCheckpoint brute_checks[] = {
            {2, 2},
            {2, 3},
            {3, 3},
        };

        for (const BruteCheckpoint& checkpoint : brute_checks) {
            if (checkpoint.m > n_ || checkpoint.n > n_) {
                continue;
            }

            const u64 brute = brute_connected_bipartite(checkpoint.m, checkpoint.n);
            const u64 got = r_value(checkpoint.m, checkpoint.n);
            if (got != brute) {
                std::cerr << "Brute-force checkpoint failed for R(" << checkpoint.m
                          << ',' << checkpoint.n << "): expected " << brute
                          << ", got " << got << '\n';
                return false;
            }

            std::cout << "Brute R(" << checkpoint.m << ',' << checkpoint.n << ") = " << got
                      << " (ok)" << '\n';
        }

        if (n_ >= 5) {
            const u64 lhs = r_value(3, 5);
            const u64 rhs = r_value(5, 3);
            if (lhs != rhs) {
                std::cerr << "Symmetry checkpoint failed: R(3,5)=" << lhs
                          << ", R(5,3)=" << rhs << '\n';
                return false;
            }
            std::cout << "Symmetry R(3,5)=R(5,3) = " << lhs << " (ok)" << '\n';
        }

        return true;
    }

private:
    int n_ = 0;
    bool allow_multithreading_ = true;
    unsigned requested_threads_ = 0;

    std::vector<std::vector<u64>> binom_;
    std::vector<u64> pow2_;
    std::vector<std::vector<u64>> r_;

    void build_binomials() {
        binom_[0][0] = 1;
        for (int n = 1; n <= n_; ++n) {
            binom_[static_cast<std::size_t>(n)][0] = 1;
            binom_[static_cast<std::size_t>(n)][static_cast<std::size_t>(n)] = 1;
            for (int k = 1; k < n; ++k) {
                binom_[static_cast<std::size_t>(n)][static_cast<std::size_t>(k)] =
                    add_mod(binom_[static_cast<std::size_t>(n - 1)]
                                  [static_cast<std::size_t>(k - 1)],
                            binom_[static_cast<std::size_t>(n - 1)]
                                  [static_cast<std::size_t>(k)]);
            }
        }
    }

    void build_powers_of_two() {
        pow2_[0] = 1;
        for (int e = 1; e <= n_ * n_; ++e) {
            pow2_[static_cast<std::size_t>(e)] =
                add_mod(pow2_[static_cast<std::size_t>(e - 1)],
                        pow2_[static_cast<std::size_t>(e - 1)]);
        }
    }

    void accumulate_previous_row_contribution(int m,
                                              int i,
                                              std::vector<u64>& out) const {
        const u64 choose_rows = binom_[static_cast<std::size_t>(m - 1)]
                                      [static_cast<std::size_t>(i - 1)];
        const int delta = m - i;

        for (int n = 1; n <= n_; ++n) {
            u64 inner = 0;
            for (int j = 0; j <= n; ++j) {
                const u64 choose_cols =
                    binom_[static_cast<std::size_t>(n)][static_cast<std::size_t>(j)];
                const u64 rij = r_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                const u64 free_part =
                    pow2_[static_cast<std::size_t>(delta * (n - j))];

                const u64 term = mul_mod(mul_mod(choose_cols, rij), free_part);
                inner = add_mod(inner, term);
            }

            out[static_cast<std::size_t>(n)] =
                add_mod(out[static_cast<std::size_t>(n)], mul_mod(choose_rows, inner));
        }
    }

    std::vector<u64> compute_base_from_previous_rows(int m) const {
        std::vector<u64> base(static_cast<std::size_t>(n_) + 1, 0);
        if (m <= 1) {
            return base;
        }

        const std::size_t workload =
            static_cast<std::size_t>(m - 1) * static_cast<std::size_t>(n_) *
            static_cast<std::size_t>(n_ + 1) / 2;
        const unsigned threads =
            choose_thread_count(allow_multithreading_, requested_threads_, workload);

        if (threads == 1) {
            for (int i = 1; i < m; ++i) {
                accumulate_previous_row_contribution(m, i, base);
            }
            return base;
        }

        std::atomic<int> next_i{1};
        std::vector<std::vector<u64>> partial(
            threads, std::vector<u64>(static_cast<std::size_t>(n_) + 1, 0));
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&, tid]() {
                while (true) {
                    const int i = next_i.fetch_add(1, std::memory_order_relaxed);
                    if (i >= m) {
                        break;
                    }
                    accumulate_previous_row_contribution(m, i, partial[tid]);
                }
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        for (unsigned tid = 0; tid < threads; ++tid) {
            for (int n = 1; n <= n_; ++n) {
                base[static_cast<std::size_t>(n)] =
                    add_mod(base[static_cast<std::size_t>(n)],
                            partial[tid][static_cast<std::size_t>(n)]);
            }
        }

        return base;
    }
};

} // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        std::cerr << "Usage: ./Euler434 [--n=100] [--single-thread] [--threads=K]" << '\n';
        return 1;
    }

    const auto t0 = std::chrono::steady_clock::now();

    Euler434Solver solver(options.n,
                          options.allow_multithreading,
                          options.requested_threads);
    solver.solve();

    if (!solver.run_checkpoints()) {
        return 1;
    }

    const u64 answer = solver.s_value(options.n);

    const auto t1 = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = t1 - t0;

    std::cout << "S(" << options.n << ") mod " << kMod << " = " << answer << '\n';
    std::cout << "Elapsed: " << elapsed.count() << " s" << '\n';

    return 0;
}
