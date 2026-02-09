#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'000ULL;

struct Options {
    int n_max = 10'000;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n-max=", options.n_max)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.n_max <= 0) {
        std::cerr << "--n-max must be positive.\n";
        return false;
    }
    return true;
}

template <typename T>
using Matrix = std::array<std::array<T, 3>, 3>;

template <typename T>
Matrix<T> zero_matrix() {
    Matrix<T> m{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m[i][j] = static_cast<T>(0);
        }
    }
    return m;
}

template <typename T>
Matrix<T> add_matrix(const Matrix<T>& a, const Matrix<T>& b, const T mod = 0) {
    Matrix<T> out = zero_matrix<T>();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (mod == 0) {
                out[i][j] = a[i][j] + b[i][j];
            } else {
                out[i][j] = (a[i][j] + b[i][j]) % mod;
            }
        }
    }
    return out;
}

template <typename T>
Matrix<T> base_matrix_for_orientation(const int from, const int to) {
    Matrix<T> m = zero_matrix<T>();
    m[from][to] = static_cast<T>(1);
    return m;
}

u64 hit_mod(const u64 k, const u64 x, const u64 y) {
    if (x == y) {
        return 0ULL;
    }
    if (x < y) {
        const u64 a = (y + kMod - x) % kMod;
        const u64 b = (x + y + kMod - 2ULL) % kMod;
        return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
    }
    const u64 a = (x + kMod - y) % kMod;
    const u64 b = (2ULL * k + kMod - x - y) % kMod;
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

u128 hit_exact(const u64 k, const u64 x, const u64 y) {
    if (x == y) {
        return 0;
    }
    if (x < y) {
        return static_cast<u128>(y - x) * static_cast<u128>(x + y - 2ULL);
    }
    return static_cast<u128>(x - y) * static_cast<u128>(2ULL * k - x - y);
}

std::string to_string_u128(u128 v) {
    if (v == 0) {
        return "0";
    }
    std::string s;
    while (v > 0) {
        s.push_back(static_cast<char>('0' + static_cast<int>(v % 10U)));
        v /= 10U;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct OrientationIndex {
    std::vector<std::tuple<int, int, int>> triples;
    std::map<std::tuple<int, int, int>, int> index_by_tuple;
};

OrientationIndex build_orientations() {
    OrientationIndex out;
    const std::vector<std::tuple<int, int, int>> all = {
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    out.triples = all;
    for (int i = 0; i < static_cast<int>(all.size()); ++i) {
        out.index_by_tuple[all[static_cast<std::size_t>(i)]] = i;
    }
    return out;
}

template <typename T>
std::vector<Matrix<T>> next_transition_counts(const std::vector<Matrix<T>>& prev,
                                              const OrientationIndex& ori,
                                              const T mod = 0) {
    std::vector<Matrix<T>> next(6, zero_matrix<T>());
    for (int id = 0; id < 6; ++id) {
        const auto [f, t, a] = ori.triples[static_cast<std::size_t>(id)];
        const int id1 = ori.index_by_tuple.at({f, a, t});
        const int id2 = ori.index_by_tuple.at({a, t, f});

        Matrix<T> cur = add_matrix(prev[static_cast<std::size_t>(id1)],
                                   prev[static_cast<std::size_t>(id2)], mod);
        if (mod == 0) {
            cur[a][f] += 1;
            cur[f][t] += 1;
            cur[t][a] += 1;
        } else {
            cur[a][f] = (cur[a][f] + 1) % mod;
            cur[f][t] = (cur[f][t] + 1) % mod;
            cur[t][a] = (cur[t][a] + 1) % mod;
        }
        next[static_cast<std::size_t>(id)] = cur;
    }
    return next;
}

u64 E_mod(const Matrix<u64>& c, const u64 k, const u64 a, const u64 b, const u64 cc) {
    std::array<u64, 3> pos{a, b, cc};
    u64 out = hit_mod(k, b, a);  // initial walk to first source rod
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (c[i][j] == 0ULL) {
                continue;
            }
            const u64 h = hit_mod(k, pos[static_cast<std::size_t>(i)],
                                  pos[static_cast<std::size_t>(j)]);
            out = (out + static_cast<u64>((static_cast<u128>(c[i][j]) * h) % kMod)) % kMod;
        }
    }
    return out;
}

u128 E_exact(const Matrix<u128>& c, const u64 k, const u64 a, const u64 b, const u64 cc) {
    std::array<u64, 3> pos{a, b, cc};
    u128 out = hit_exact(k, b, a);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (c[i][j] == 0U) {
                continue;
            }
            out += c[i][j] * hit_exact(k, pos[static_cast<std::size_t>(i)],
                                       pos[static_cast<std::size_t>(j)]);
        }
    }
    return out;
}

u64 solve_sum_last9(const int n_max) {
    const OrientationIndex ori = build_orientations();
    const int canonical = ori.index_by_tuple.at({0, 2, 1});  // move from rod 0 to rod 2 using rod 1

    std::vector<Matrix<u64>> cur_mod(6, zero_matrix<u64>());
    for (int id = 0; id < 6; ++id) {
        const auto [f, t, a] = ori.triples[static_cast<std::size_t>(id)];
        (void)a;
        cur_mod[static_cast<std::size_t>(id)] = base_matrix_for_orientation<u64>(f, t);
    }

    u64 p3 = 1ULL;
    u64 p6 = 1ULL;
    u64 p9 = 1ULL;
    u64 p10 = 1ULL;

    u64 total = 0ULL;
    for (int n = 1; n <= n_max; ++n) {
        if (n > 1) {
            cur_mod = next_transition_counts(cur_mod, ori, kMod);
        }

        p3 = static_cast<u64>((static_cast<u128>(p3) * 3ULL) % kMod);
        p6 = static_cast<u64>((static_cast<u128>(p6) * 6ULL) % kMod);
        p9 = static_cast<u64>((static_cast<u128>(p9) * 9ULL) % kMod);
        p10 = static_cast<u64>((static_cast<u128>(p10) * 10ULL) % kMod);

        const u64 e = E_mod(cur_mod[static_cast<std::size_t>(canonical)], p10, p3, p6, p9);
        total += e;
        total %= kMod;
    }
    return total;
}

bool run_checkpoints() {
    const OrientationIndex ori = build_orientations();
    const int canonical = ori.index_by_tuple.at({0, 2, 1});

    std::vector<Matrix<u128>> cur(6, zero_matrix<u128>());
    for (int id = 0; id < 6; ++id) {
        const auto [f, t, a] = ori.triples[static_cast<std::size_t>(id)];
        (void)a;
        cur[static_cast<std::size_t>(id)] = base_matrix_for_orientation<u128>(f, t);
    }

    // n = 2
    cur = next_transition_counts(cur, ori, static_cast<u128>(0));
    if (E_exact(cur[static_cast<std::size_t>(canonical)], 5ULL, 1ULL, 3ULL, 5ULL) !=
        static_cast<u128>(60ULL)) {
        std::cerr << "Checkpoint failed: E(2,5,1,3,5)\n";
        return false;
    }

    // n = 3
    cur = next_transition_counts(cur, ori, static_cast<u128>(0));
    if (E_exact(cur[static_cast<std::size_t>(canonical)], 20ULL, 4ULL, 9ULL, 17ULL) !=
        static_cast<u128>(2358ULL)) {
        std::cerr << "Checkpoint failed: E(3,20,4,9,17)\n";
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
        return 1;
    }

    const u64 ans = solve_sum_last9(options.n_max);
    std::cout << std::setw(9) << std::setfill('0') << ans << '\n';
    return 0;
}
