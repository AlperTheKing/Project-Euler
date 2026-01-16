#include <algorithm>
#include <cstdint>
#include <functional>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = uint32_t;
using u64 = uint64_t;

constexpr u32 kMod = 1000000000u;
constexpr int kDefaultM = 10;
constexpr long long kDefaultN = 1000000000000LL;

struct Matrix {
    int n = 0;
    std::vector<u32> a;

    Matrix() = default;
    explicit Matrix(int n_) : n(n_), a(static_cast<size_t>(n_) * n_, 0) {}

    u32 &operator()(int r, int c) { return a[static_cast<size_t>(r) * n + c]; }
    const u32 &operator()(int r, int c) const { return a[static_cast<size_t>(r) * n + c]; }
};

Matrix identity(int n) {
    Matrix m(n);
    for (int i = 0; i < n; ++i) m(i, i) = 1;
    return m;
}

Matrix transpose(const Matrix &m) {
    Matrix t(m.n);
    for (int i = 0; i < m.n; ++i) {
        for (int j = 0; j < m.n; ++j) {
            t(j, i) = m(i, j);
        }
    }
    return t;
}

Matrix add(const Matrix &a, const Matrix &b) {
    Matrix c(a.n);
    const size_t size = a.a.size();
    for (size_t i = 0; i < size; ++i) {
        u32 val = a.a[i] + b.a[i];
        if (val >= kMod) val -= kMod;
        c.a[i] = val;
    }
    return c;
}

Matrix subtract(const Matrix &a, const Matrix &b) {
    Matrix c(a.n);
    const size_t size = a.a.size();
    for (size_t i = 0; i < size; ++i) {
        u32 val = (a.a[i] >= b.a[i]) ? (a.a[i] - b.a[i]) : (a.a[i] + kMod - b.a[i]);
        c.a[i] = val;
    }
    return c;
}

Matrix multiply(const Matrix &a, const Matrix &b, int threads) {
    const int n = a.n;
    Matrix c(n);
    Matrix bt = transpose(b);

    auto worker = [&](int row_start, int row_end) {
        for (int i = row_start; i < row_end; ++i) {
            const u32 *row_a = &a.a[static_cast<size_t>(i) * n];
            for (int j = 0; j < n; ++j) {
                const u32 *row_b = &bt.a[static_cast<size_t>(j) * n];
                unsigned __int128 sum = 0;
                for (int k = 0; k < n; ++k) {
                    sum += static_cast<unsigned __int128>(row_a[k]) * row_b[k];
                }
                c(i, j) = static_cast<u32>(sum % kMod);
            }
        }
    };

    if (threads <= 1 || n < 64) {
        worker(0, n);
        return c;
    }

    int use_threads = std::min(threads, n);
    int chunk = (n + use_threads - 1) / use_threads;
    std::vector<std::thread> pool;
    pool.reserve(use_threads);
    for (int t = 0; t < use_threads; ++t) {
        int start = t * chunk;
        int end = std::min(n, start + chunk);
        if (start >= end) break;
        pool.emplace_back(worker, start, end);
    }
    for (auto &th : pool) th.join();
    return c;
}

struct Triple {
    int a = 0;
    int b = 0;
    int c = 0;
};

struct Basis {
    int p = 0;
    int dim = 0;
    int idx_c_pow = 0;
    std::vector<Triple> triples;
    std::vector<std::vector<int>> idx;
};

Basis build_basis(int p) {
    Basis basis;
    basis.p = p;
    basis.idx.assign(p + 1, std::vector<int>(p + 1, -1));
    for (int a = 0; a <= p; ++a) {
        for (int b = 0; b <= p - a; ++b) {
            int c = p - a - b;
            basis.idx[a][b] = static_cast<int>(basis.triples.size());
            basis.triples.push_back({a, b, c});
        }
    }
    basis.dim = static_cast<int>(basis.triples.size());
    basis.idx_c_pow = basis.idx[0][0];
    return basis;
}

std::vector<std::vector<u64>> build_combinations(int p) {
    std::vector<std::vector<u64>> comb(p + 1, std::vector<u64>(p + 1, 0));
    comb[0][0] = 1;
    for (int n = 1; n <= p; ++n) {
        comb[n][0] = comb[n][n] = 1;
        for (int k = 1; k < n; ++k) {
            comb[n][k] = comb[n - 1][k - 1] + comb[n - 1][k];
        }
    }
    return comb;
}

Matrix build_allowed(const Basis &basis, const std::vector<std::vector<u64>> &comb) {
    Matrix allowed(basis.dim);
    for (int idx_new = 0; idx_new < basis.dim; ++idx_new) {
        const Triple &t = basis.triples[idx_new];
        int i = t.a;
        int j = t.b;
        int k = t.c;
        for (int u = 0; u <= k; ++u) {
            for (int v = 0; v <= k - u; ++v) {
                int w = k - u - v;
                int a_old = u;
                int b_old = v + i;
                int idx_old = basis.idx[a_old][b_old];
                u64 coeff = comb[k][u] * comb[k - u][v];
                u32 add = static_cast<u32>(coeff % kMod);
                u32 &cell = allowed(idx_new, idx_old);
                cell = static_cast<u32>((cell + add) % kMod);
            }
        }
    }
    return allowed;
}

Matrix build_forbidden(const Basis &basis) {
    Matrix forbidden(basis.dim);
    for (int idx_new = 0; idx_new < basis.dim; ++idx_new) {
        const Triple &t = basis.triples[idx_new];
        if (t.c != 0) continue;
        int idx_old = basis.idx[0][t.a];
        forbidden(idx_new, idx_old) = 1;
    }
    return forbidden;
}

struct PowResult {
    Matrix value;
    Matrix deriv;
};

PowResult pow_with_derivative(const Matrix &base, const Matrix &base_deriv, u64 exp, int threads) {
    const int n = base.n;
    Matrix result = identity(n);
    Matrix result_deriv(n);
    Matrix cur = base;
    Matrix cur_deriv = base_deriv;

    u64 e = exp;
    while (e > 0) {
        if (e & 1ULL) {
            Matrix tmp1 = multiply(result_deriv, cur, threads);
            Matrix tmp2 = multiply(result, cur_deriv, threads);
            Matrix next_deriv = add(tmp1, tmp2);
            Matrix next = multiply(result, cur, threads);
            result = std::move(next);
            result_deriv = std::move(next_deriv);
        }
        e >>= 1ULL;
        if (e == 0) break;
        Matrix tmp1 = multiply(cur_deriv, cur, threads);
        Matrix tmp2 = multiply(cur, cur_deriv, threads);
        Matrix next_deriv = add(tmp1, tmp2);
        Matrix next = multiply(cur, cur, threads);
        cur = std::move(next);
        cur_deriv = std::move(next_deriv);
    }
    return {std::move(result), std::move(result_deriv)};
}

struct Solver {
    int m = 0;
    int p = 0;
    Basis basis;
    Matrix allowed;
    Matrix forbidden;
    Matrix combined;

    explicit Solver(int m_) : m(m_), p(2 * m_), basis(build_basis(2 * m_)) {
        auto comb = build_combinations(p);
        allowed = build_allowed(basis, comb);
        forbidden = build_forbidden(basis);
        combined = subtract(allowed, forbidden);  // y = -1
    }

    u32 solve(long long n, int threads) const {
        if (n <= 1) return 1;
        u64 internal = (n >= 2) ? static_cast<u64>(n - 2) : 0;
        PowResult pow = (internal == 0)
                            ? PowResult{identity(basis.dim), Matrix(basis.dim)}
                            : pow_with_derivative(combined, forbidden, internal, threads);

        const int idx0 = basis.idx_c_pow;
        std::vector<u32> col(pow.value.n);
        std::vector<u32> col_deriv(pow.value.n);
        for (int i = 0; i < pow.value.n; ++i) {
            col[i] = pow.value(i, idx0);
            col_deriv[i] = pow.deriv(i, idx0);
        }

        unsigned __int128 sum_val = 0;
        unsigned __int128 sum_deriv = 0;
        const u32 *row = &allowed.a[static_cast<size_t>(idx0) * allowed.n];
        for (int j = 0; j < allowed.n; ++j) {
            sum_val += static_cast<unsigned __int128>(row[j]) * col[j];
            sum_deriv += static_cast<unsigned __int128>(row[j]) * col_deriv[j];
        }
        u32 f_val = static_cast<u32>(sum_val % kMod);
        u32 f_deriv = static_cast<u32>(sum_deriv % kMod);
        u32 ans = f_val + f_deriv;
        if (ans >= kMod) ans -= kMod;
        return ans;
    }
};

std::vector<int> build_path_masks(int n) {
    std::vector<int> masks;
    if (n <= 0) return masks;

    std::function<void(int, int)> dfs = [&](int pos, int mask) {
        if (pos == n) {
            masks.push_back(mask);
            return;
        }
        for (int step = 1; step <= 3; ++step) {
            int next = pos + step;
            if (next > n) continue;
            int next_mask = mask | (1 << (next - 1));
            dfs(next, next_mask);
        }
    };

    dfs(1, 1);
    return masks;
}

u64 brute_count(int m, int n) {
    if (n <= 1) return 1;
    std::vector<int> masks = build_path_masks(n);
    int legs = 2 * m;
    int max_mask = 1 << n;
    std::vector<u64> dp(max_mask, 0);
    dp[0] = 1;
    for (int step = 0; step < legs; ++step) {
        std::vector<u64> next(max_mask, 0);
        for (int mask = 0; mask < max_mask; ++mask) {
            u64 val = dp[mask];
            if (val == 0) continue;
            for (int pmask : masks) {
                next[mask | pmask] += val;
            }
        }
        dp.swap(next);
    }

    int internal = std::max(0, n - 2);
    int internal_mask = ((1 << n) - 1) ^ 1 ^ (1 << (n - 1));
    u64 total = 0;
    for (int mask = 0; mask < max_mask; ++mask) {
        u64 val = dp[mask];
        if (val == 0) continue;
        int visited_internal = __builtin_popcount(mask & internal_mask);
        int unvisited = internal - visited_internal;
        if (unvisited <= 1) total += val;
    }
    return total % kMod;
}

void run_validation() {
    struct Check {
        int m;
        long long n;
        u32 expected;
    };
    const Check checks[] = {
        {1, 3, 4},
        {1, 4, 15},
        {1, 5, 46},
        {2, 3, 16},
        {2, 100, 429619151},
    };
    for (const auto &chk : checks) {
        Solver solver(chk.m);
        u32 got = solver.solve(chk.n, 1);
        if (got != chk.expected) {
            std::cerr << "Validation failed for m=" << chk.m << ", n=" << chk.n
                      << ": got " << got << ", expected " << chk.expected << '\n';
            std::exit(1);
        }
    }

    for (int m = 1; m <= 2; ++m) {
        Solver solver(m);
        for (int n = 3; n <= 7; ++n) {
            u32 fast = solver.solve(n, 1);
            u32 brute = static_cast<u32>(brute_count(m, n));
            if (fast != brute) {
                std::cerr << "Brute validation failed for m=" << m << ", n=" << n
                          << ": got " << fast << ", expected " << brute << '\n';
                std::exit(1);
            }
        }
    }
    std::cerr << "Validation checkpoints passed.\n";
}

void usage(const char *argv0) {
    std::cerr << "Usage: " << argv0 << " [-m M] [-n N] [-t THREADS] [--no-validate]\n";
}

}  // namespace

int main(int argc, char **argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int m = kDefaultM;
    long long n = kDefaultN;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;
    bool validate = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-m" && i + 1 < argc) {
            m = std::stoi(argv[++i]);
        } else if (arg == "-n" && i + 1 < argc) {
            n = std::stoll(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            threads = std::max(1, std::stoi(argv[++i]));
        } else if (arg == "--no-validate") {
            validate = false;
        } else if (arg == "--validate") {
            validate = true;
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    if (validate) run_validation();

    Solver solver(m);
    u32 answer = solver.solve(n, threads);
    std::cout << answer << '\n';
    return 0;
}
