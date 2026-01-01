#include <algorithm>
#include <cstdint>
#include <iostream>
#include <functional>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr int MOD = 998244353;

int mod_pow(int base, int exp) {
    long long res = 1;
    long long cur = base;
    while (exp > 0) {
        if (exp & 1) res = (res * cur) % MOD;
        cur = (cur * cur) % MOD;
        exp >>= 1;
    }
    return static_cast<int>(res);
}

int add_mod(int a, int b) {
    int v = a + b;
    if (v >= MOD) v -= MOD;
    return v;
}

int sub_mod(int a, int b) {
    int v = a - b;
    if (v < 0) v += MOD;
    return v;
}

int compute_c(int n, uint64_t b) {
    if (n == 1) return static_cast<int>((b + 1) % MOD);
    if (b == 0) return 0;

    vector<int> positions;
    for (int i = 63; i >= 0; --i) {
        if ((b >> i) & 1ULL) positions.push_back(i);
    }
    int m = static_cast<int>(positions.size());
    if (m == 0) return 0;

    vector<int> seg(m, 0);
    for (int i = 0; i + 1 < m; ++i) {
        seg[i] = positions[i] - positions[i + 1] - 1;
    }
    seg[m - 1] = positions.back();

    vector<int> dims(m, 1);
    vector<int> stride(m, 1);
    int P = 1;
    for (int i = 0; i < m; ++i) {
        dims[i] = seg[i] + 1;
        if (i > 0) stride[i] = stride[i - 1] * dims[i - 1];
        P *= dims[i];
    }

    int max_seg = 0;
    for (int s : seg) max_seg = max(max_seg, s);
    vector<int> fact(max_seg + 1, 1), inv_fact(max_seg + 1, 1);
    for (int i = 1; i <= max_seg; ++i) {
        fact[i] = static_cast<int>((1LL * fact[i - 1] * i) % MOD);
    }
    if (max_seg > 0) {
        inv_fact[max_seg] = mod_pow(fact[max_seg], MOD - 2);
        for (int i = max_seg; i > 0; --i) {
            inv_fact[i - 1] = static_cast<int>((1LL * inv_fact[i] * i) % MOD);
        }
    }
    auto comb = [&](int nn, int kk) -> int {
        if (kk < 0 || kk > nn) return 0;
        return static_cast<int>(
            (1LL * fact[nn] * inv_fact[kk] % MOD) * inv_fact[nn - kk] % MOD);
    };

    vector<vector<int>> combs(m);
    for (int i = 0; i < m; ++i) {
        combs[i].assign(seg[i] + 1, 0);
        for (int k = 0; k <= seg[i]; ++k) combs[i][k] = comb(seg[i], k);
    }

    vector<int> size_c(P, 0);
    vector<uint32_t> nz_mask(P, 0);
    vector<int> counts(m, 0);
    function<void(int, int)> enumerate_counts = [&](int idx, int dim) {
        if (dim == m) {
            long long prod = 1;
            uint32_t mask = 0;
            for (int i = 0; i < m; ++i) {
                int c = counts[i];
                prod = (prod * combs[i][c]) % MOD;
                if (c > 0) mask |= (1u << i);
            }
            size_c[idx] = static_cast<int>(prod);
            nz_mask[idx] = mask;
            return;
        }
        int step = stride[dim];
        for (int c = 0; c <= seg[dim]; ++c) {
            counts[dim] = c;
            enumerate_counts(idx + c * step, dim + 1);
        }
    };
    enumerate_counts(0, 0);

    int S = 1 << m;
    vector<int> boundary(S, m);
    for (int s = 0; s < S; ++s) {
        int t = m;
        for (int i = 0; i < m; ++i) {
            if ((s & (1 << i)) == 0) {
                t = i;
                break;
            }
        }
        boundary[s] = t;
    }
    vector<uint32_t> prefix_mask(m + 1, 0);
    for (int t = 1; t <= m; ++t) prefix_mask[t] = prefix_mask[t - 1] | (1u << (t - 1));

    vector<uint32_t> prefix_by_s(S, 0);
    vector<int> comp_s(S, 0);
    int full_mask = S - 1;
    for (int s = 0; s < S; ++s) {
        prefix_by_s[s] = prefix_mask[boundary[s]];
        comp_s[s] = full_mask ^ s;
    }

    vector<vector<int>> E(m);
    int max_dim = 1;
    for (int i = 0; i < m; ++i) {
        int di = dims[i];
        max_dim = max(max_dim, di);
        E[i].assign(di * di, 0);
        for (int a = 0; a < di; ++a) {
            for (int b = 0; b < di; ++b) {
                if (seg[i] - a >= b) E[i][a * di + b] = comb(seg[i] - a, b);
            }
        }
    }

    size_t total_size = static_cast<size_t>(S) * static_cast<size_t>(P);
    vector<int> f(total_size, 0);
    vector<int> g(total_size, 0);

    unsigned hw = thread::hardware_concurrency();
    int threads = hw == 0 ? 1 : static_cast<int>(hw);
    if (total_size < 1'000'000) threads = 1;

    auto parallel_for = [&](int start, int end, int pieces, const auto& fn) {
        if (pieces <= 1 || end - start <= 1) {
            fn(0, start, end);
            return;
        }
        int total = end - start;
        int chunk = (total + pieces - 1) / pieces;
        vector<thread> workers;
        for (int t = 0; t < pieces; ++t) {
            int s = start + t * chunk;
            int e = min(end, s + chunk);
            if (s >= e) break;
            workers.emplace_back([=, &fn]() { fn(t, s, e); });
        }
        for (auto& th : workers) th.join();
    };

    parallel_for(0, S, threads, [&](int, int s_begin, int s_end) {
        for (int s = s_begin; s < s_end; ++s) {
            uint32_t pref = prefix_by_s[s];
            size_t base = static_cast<size_t>(s) * P;
            for (int c = 0; c < P; ++c) {
                if ((nz_mask[c] & pref) != 0u) continue;
                if (s == 0 && c == 0) continue;
                f[base + c] = 1;
            }
        }
    });

    auto zeta_transform = [&](vector<int>& data) {
        for (int bit = 0; bit < m; ++bit) {
            int step = 1 << bit;
            int block = step << 1;
            parallel_for(0, P, threads, [&](int, int c_begin, int c_end) {
                for (int s0 = 0; s0 < S; s0 += block) {
                    for (int k = 0; k < step; ++k) {
                        size_t idx1 = static_cast<size_t>(s0 + k) * P + c_begin;
                        size_t idx2 = static_cast<size_t>(s0 + k + step) * P + c_begin;
                        for (int c = c_begin; c < c_end; ++c) {
                            int v = data[idx2 + (c - c_begin)] + data[idx1 + (c - c_begin)];
                            if (v >= MOD) v -= MOD;
                            data[idx2 + (c - c_begin)] = v;
                        }
                    }
                }
            });
        }
    };

    auto apply_E_transform = [&](vector<int>& data) {
        parallel_for(0, S, threads, [&](int, int s_begin, int s_end) {
            vector<int> tmp(max_dim, 0);
            for (int s = s_begin; s < s_end; ++s) {
                size_t base = static_cast<size_t>(s) * P;
                for (int i = 0; i < m; ++i) {
                    int di = dims[i];
                    if (di == 1) continue;
                    int step = stride[i];
                    int block = step * di;
                    const int* mat = E[i].data();
                    for (int start = 0; start < P; start += block) {
                        for (int off = 0; off < step; ++off) {
                            for (int k = 0; k < di; ++k) {
                                tmp[k] = data[base + start + off + k * step];
                            }
                            for (int a = 0; a < di; ++a) {
                                long long sum = 0;
                                const int* row = mat + a * di;
                                for (int b = 0; b < di; ++b) {
                                    sum += 1LL * row[b] * tmp[b];
                                }
                                data[base + start + off + a * step] = static_cast<int>(sum % MOD);
                            }
                        }
                    }
                }
            }
        });
    };

    auto total_weight = [&](const vector<int>& data) -> int {
        vector<long long> partial(threads, 0);
        parallel_for(0, S, threads, [&](int tid, int s_begin, int s_end) {
            long long sum = 0;
            for (int s = s_begin; s < s_end; ++s) {
                size_t base = static_cast<size_t>(s) * P;
                for (int c = 0; c < P; ++c) {
                    if (data[base + c] == 0) continue;
                    sum += 1LL * data[base + c] * size_c[c];
                    if (sum >= (1LL << 62)) sum %= MOD;
                }
            }
            partial[tid] = sum % MOD;
        });
        long long sum = 0;
        for (long long v : partial) {
            sum += v;
            if (sum >= (1LL << 62)) sum %= MOD;
        }
        return static_cast<int>(sum % MOD);
    };

    auto build_next = [&](int total, const vector<int>& src, vector<int>& dst) {
        parallel_for(0, S, threads, [&](int, int s_begin, int s_end) {
            for (int s = s_begin; s < s_end; ++s) {
                size_t base = static_cast<size_t>(s) * P;
                size_t base_comp = static_cast<size_t>(comp_s[s]) * P;
                uint32_t pref = prefix_by_s[s];
                for (int c = 0; c < P; ++c) {
                    if ((nz_mask[c] & pref) != 0u || (s == 0 && c == 0)) {
                        dst[base + c] = 0;
                        continue;
                    }
                    dst[base + c] = sub_mod(total, src[base_comp + c]);
                }
            }
        });
    };

    for (int step = 0; step < n - 1; ++step) {
        int total = total_weight(f);
        zeta_transform(f);
        apply_E_transform(f);
        build_next(total, f, g);
        f.swap(g);
    }

    return total_weight(f);
}

}  // namespace

int main() {
    struct Check {
        int n;
        uint64_t b;
        int expected;
    };

    const Check checks[] = {
        {3, 4, 18},
        {10, 6, 2496120},
        {100, 200, 268159379},
    };

    for (const auto& chk : checks) {
        int got = compute_c(chk.n, chk.b);
        if (got != chk.expected) {
            cerr << "Validation failure: c(" << chk.n << ", " << chk.b
                 << ") = " << got << ", expected " << chk.expected << '\n';
            return 1;
        }
    }

    cout << compute_c(123, 123456789ULL) << '\n';
    return 0;
}
