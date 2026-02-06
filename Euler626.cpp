#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

constexpr std::uint64_t MOD = 1001001011ULL;
constexpr int MAX_N = 20;
constexpr int MAX_VARS = 2 * MAX_N;

struct CycleType {
    std::array<int, MAX_N + 1> mult{};
    std::vector<int> valuations;  // v2 of each cycle length
    int cycles = 0;
    cpp_int perm_count = 0;       // number of permutations with this cycle type
};

int v2(int x) {
    int c = 0;
    while ((x & 1) == 0) {
        x >>= 1;
        ++c;
    }
    return c;
}

std::vector<cpp_int> factorials_up_to(int n) {
    std::vector<cpp_int> fact(static_cast<std::size_t>(n + 1), 1);
    for (int i = 1; i <= n; ++i) fact[static_cast<std::size_t>(i)] = fact[static_cast<std::size_t>(i - 1)] * i;
    return fact;
}

cpp_int ipow_int(int base, int exp) {
    cpp_int result = 1;
    cpp_int b = base;
    int e = exp;
    while (e > 0) {
        if (e & 1) result *= b;
        b *= b;
        e >>= 1;
    }
    return result;
}

void build_cycle_types_rec(int n,
                           int remaining,
                           int min_part,
                           std::vector<int>& current,
                           const std::vector<cpp_int>& fact,
                           std::vector<CycleType>& out) {
    if (remaining == 0) {
        CycleType t;
        t.cycles = static_cast<int>(current.size());
        for (int len : current) {
            ++t.mult[static_cast<std::size_t>(len)];
            t.valuations.push_back(v2(len));
        }

        cpp_int denom = 1;
        for (int len = 1; len <= n; ++len) {
            int m = t.mult[static_cast<std::size_t>(len)];
            if (m == 0) continue;
            denom *= ipow_int(len, m);
            denom *= fact[static_cast<std::size_t>(m)];
        }
        t.perm_count = fact[static_cast<std::size_t>(n)] / denom;
        out.push_back(std::move(t));
        return;
    }

    for (int part = min_part; part <= remaining; ++part) {
        current.push_back(part);
        build_cycle_types_rec(n, remaining - part, part, current, fact, out);
        current.pop_back();
    }
}

std::vector<CycleType> build_cycle_types(int n, const std::vector<cpp_int>& fact) {
    std::vector<CycleType> types;
    std::vector<int> cur;
    build_cycle_types_rec(n, n, 1, cur, fact, types);
    return types;
}

struct DSU {
    std::array<int, MAX_VARS> parent{};
    std::array<unsigned char, MAX_VARS> rank{};
    std::array<unsigned char, MAX_VARS> forced{};

    void init(int n) {
        for (int i = 0; i < n; ++i) {
            parent[static_cast<std::size_t>(i)] = i;
            rank[static_cast<std::size_t>(i)] = 0;
            forced[static_cast<std::size_t>(i)] = 0;
        }
    }

    int find(int x) {
        int r = x;
        while (parent[static_cast<std::size_t>(r)] != r) r = parent[static_cast<std::size_t>(r)];
        while (parent[static_cast<std::size_t>(x)] != x) {
            int p = parent[static_cast<std::size_t>(x)];
            parent[static_cast<std::size_t>(x)] = r;
            x = p;
        }
        return r;
    }

    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (rank[static_cast<std::size_t>(a)] < rank[static_cast<std::size_t>(b)]) std::swap(a, b);
        parent[static_cast<std::size_t>(b)] = a;
        forced[static_cast<std::size_t>(a)] |= forced[static_cast<std::size_t>(b)];
        if (rank[static_cast<std::size_t>(a)] == rank[static_cast<std::size_t>(b)]) ++rank[static_cast<std::size_t>(a)];
    }

    void set_forced(int a) {
        int r = find(a);
        forced[static_cast<std::size_t>(r)] = 1;
    }
};

int parity_rank(const std::vector<int>& rows_v2, const std::vector<int>& cols_v2) {
    const int r = static_cast<int>(rows_v2.size());
    const int s = static_cast<int>(cols_v2.size());
    const int vars = r + s;
    DSU dsu;
    dsu.init(vars);

    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < s; ++j) {
            if (rows_v2[static_cast<std::size_t>(i)] == cols_v2[static_cast<std::size_t>(j)]) {
                dsu.unite(i, r + j);                 // A_i = B_j
            } else if (rows_v2[static_cast<std::size_t>(i)] > cols_v2[static_cast<std::size_t>(j)]) {
                dsu.set_forced(i);                   // A_i = 0
            } else {
                dsu.set_forced(r + j);               // B_j = 0
            }
        }
    }

    std::array<unsigned char, MAX_VARS> seen{};
    int free_components = 0;
    for (int v = 0; v < vars; ++v) {
        int root = dsu.find(v);
        if (seen[static_cast<std::size_t>(root)]) continue;
        seen[static_cast<std::size_t>(root)] = 1;
        if (!dsu.forced[static_cast<std::size_t>(root)]) ++free_components;
    }
    return vars - free_components;
}

int orbit_count(const CycleType& a, const CycleType& b, const int gcd_tbl[MAX_N + 1][MAX_N + 1], int n) {
    int o = 0;
    for (int p = 1; p <= n; ++p) {
        int mp = a.mult[static_cast<std::size_t>(p)];
        if (mp == 0) continue;
        for (int q = 1; q <= n; ++q) {
            int mq = b.mult[static_cast<std::size_t>(q)];
            if (mq == 0) continue;
            o += mp * mq * gcd_tbl[p][q];
        }
    }
    return o;
}

cpp_int compute_c_exact(int n, unsigned threads) {
    std::vector<cpp_int> fact = factorials_up_to(n);
    std::vector<CycleType> types = build_cycle_types(n, fact);

    int gcd_tbl[MAX_N + 1][MAX_N + 1]{};
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) gcd_tbl[i][j] = std::gcd(i, j);
    }

    const std::size_t tcount = types.size();
    if (threads == 0) threads = 1;
    threads = std::min<unsigned>(threads, static_cast<unsigned>(tcount));
    if (threads == 0) threads = 1;

    std::atomic<std::size_t> next_i{0};
    std::vector<cpp_int> partial(static_cast<std::size_t>(threads), 0);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    for (unsigned tid = 0; tid < threads; ++tid) {
        pool.emplace_back([&, tid]() {
            cpp_int local = 0;
            while (true) {
                std::size_t i = next_i.fetch_add(1, std::memory_order_relaxed);
                if (i >= tcount) break;
                const CycleType& a = types[i];
                for (std::size_t j = 0; j < tcount; ++j) {
                    const CycleType& b = types[j];
                    int o = orbit_count(a, b, gcd_tbl, n);
                    int rk = parity_rank(a.valuations, b.valuations);
                    int exp = o - rk;
                    assert(exp >= 0);
                    local += a.perm_count * b.perm_count * (cpp_int(1) << exp);
                }
            }
            partial[static_cast<std::size_t>(tid)] = std::move(local);
        });
    }

    for (auto& th : pool) th.join();

    cpp_int total = 0;
    for (const auto& p : partial) total += p;

    cpp_int denom = fact[static_cast<std::size_t>(n)] * fact[static_cast<std::size_t>(n)];
    assert(total % denom == 0);
    return total / denom;
}

std::uint64_t to_u64(const cpp_int& x) {
    return x.convert_to<std::uint64_t>();
}

void run_checks(unsigned threads) {
    assert(to_u64(compute_c_exact(1, 1)) == 1ULL);
    assert(to_u64(compute_c_exact(3, 1)) == 3ULL);
    assert(to_u64(compute_c_exact(5, threads)) == 39ULL);
    assert(to_u64(compute_c_exact(8, threads)) == 656108ULL);
}

}  // namespace

int main() {
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;

    run_checks(threads);

    cpp_int c20 = compute_c_exact(20, threads);
    std::uint64_t answer = to_u64(c20 % MOD);
    std::cout << answer << "\n";
    return 0;
}
