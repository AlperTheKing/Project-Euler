#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <queue>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

constexpr int MOD = 1'000'000'007;

long long mod_pow(long long a, long long e) {
    long long r = 1 % MOD;
    a %= MOD;
    while (e > 0) {
        if (e & 1LL) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1LL;
    }
    return r;
}

struct Block {
    int len;
    int tiles;
    int count;
};

struct DFA {
    std::vector<std::array<int, 5>> trans;
    std::vector<int> accept0;
    std::vector<int> accept1;
};

DFA build_dfa() {
    auto id = [](int a, int b, int p) { return p * 25 + a * 5 + b; };
    const int nfa_states = 50;

    std::array<std::array<uint64_t, 5>, nfa_states> nfa{};

    for (int p = 0; p <= 1; ++p) {
        for (int a = 0; a <= 4; ++a) {
            for (int b = 0; b <= 4; ++b) {
                const int sid = id(a, b, p);
                for (int c = 0; c <= 4; ++c) {
                    uint64_t mask = 0;
                    if (c >= a + b) {
                        const int remaining = c - a - b;
                        for (int use_pair = 0; use_pair <= 1; ++use_pair) {
                            if (use_pair && p) continue;
                            if (use_pair && remaining < 2) continue;
                            const int rem2 = remaining - 2 * use_pair;
                            for (int q = 0; q <= 4; ++q) {
                                if (q > rem2) break;
                                if ((rem2 - q) % 3 != 0) continue;
                                const int ns = id(q, a, p || use_pair);
                                mask |= 1ULL << ns;
                            }
                        }
                    }
                    nfa[sid][c] = mask;
                }
            }
        }
    }

    const uint64_t start_mask = 1ULL << id(0, 0, 0);
    std::unordered_map<uint64_t, int> index;
    std::vector<uint64_t> masks;
    std::queue<uint64_t> q;
    index[start_mask] = 0;
    masks.push_back(start_mask);
    q.push(start_mask);

    std::vector<std::array<int, 5>> trans;

    while (!q.empty()) {
        const uint64_t mask = q.front();
        q.pop();
        const int cur = index[mask];
        if (cur >= static_cast<int>(trans.size())) {
            trans.push_back({0, 0, 0, 0, 0});
        }
        for (int c = 0; c <= 4; ++c) {
            uint64_t next_mask = 0;
            uint64_t tmp = mask;
            while (tmp) {
                const int i = __builtin_ctzll(tmp);
                tmp &= tmp - 1;
                next_mask |= nfa[i][c];
            }
            auto it = index.find(next_mask);
            int nxt;
            if (it == index.end()) {
                nxt = static_cast<int>(index.size());
                index[next_mask] = nxt;
                masks.push_back(next_mask);
                q.push(next_mask);
            } else {
                nxt = it->second;
            }
            trans[cur][c] = nxt;
        }
    }

    const uint64_t acc0_bit = 1ULL << id(0, 0, 0);
    const uint64_t acc1_bit = 1ULL << id(0, 0, 1);

    std::vector<int> accept0;
    std::vector<int> accept1;
    for (int i = 0; i < static_cast<int>(masks.size()); ++i) {
        if (masks[i] & acc0_bit) accept0.push_back(i);
        if (masks[i] & acc1_bit) accept1.push_back(i);
    }

    DFA dfa;
    dfa.trans = std::move(trans);
    dfa.accept0 = std::move(accept0);
    dfa.accept1 = std::move(accept1);
    return dfa;
}

class MahjongCounter {
public:
    MahjongCounter(int max_t, unsigned threads)
        : max_t_(max_t),
          max_tiles_(3 * max_t + 2),
          max_len_(3 * max_t + 2),
          max_blocks_(max_t + 1),
          stride_(max_tiles_ + 1),
          threads_(std::max(1u, threads)),
          dfa_(build_dfa()) {
        build_factorials();
        build_block_counts();
        build_block_dp();
    }

    long long solve(long long n, long long s, int t) const {
        if (s == 0) return 0;
        if (t > max_t_) return 0;

        const int tiles_max = 3 * t + 2;
        const int len_max = tiles_max;
        const int blocks_max = t + 1;

        std::vector<int> A(tiles_max + 1, 0);
        std::vector<int> B(tiles_max + 1, 0);

        for (int k = 0; k <= blocks_max; ++k) {
            for (int L = 0; L <= len_max; ++L) {
                const long long N = n - L + 1;
                if (N < k) continue;
                const int comb = comb_small(N, k);
                if (comb == 0) continue;
                const int base = L * stride_;
                const std::vector<int>& dp0 = dp0_[k];
                const std::vector<int>& dp1 = dp1_[k];
                for (int M = 0; M <= tiles_max; ++M) {
                    const int v0 = dp0[base + M];
                    if (v0) {
                        A[M] = (A[M] + static_cast<long long>(v0) * comb) % MOD;
                    }
                    const int v1 = dp1[base + M];
                    if (v1) {
                        B[M] = (B[M] + static_cast<long long>(v1) * comb) % MOD;
                    }
                }
            }
        }

        std::vector<int> G(t + 1, 0);
        std::vector<int> H(t + 1, 0);
        for (int k = 0; k <= t; ++k) {
            const int idxG = 3 * k;
            const int idxH = 3 * k + 2;
            if (idxG <= tiles_max) G[k] = A[idxG];
            if (idxH <= tiles_max) H[k] = B[idxH];
        }

        std::vector<int> Gpow = poly_pow(G, s - 1, t);
        long long total = 0;
        for (int i = 0; i <= t; ++i) {
            total = (total + static_cast<long long>(H[i]) * Gpow[t - i]) % MOD;
        }
        total = total * (s % MOD) % MOD;
        return total;
    }

private:
    int max_t_;
    int max_tiles_;
    int max_len_;
    int max_blocks_;
    int stride_;
    unsigned threads_;
    DFA dfa_;

    std::vector<int> fact_;
    std::vector<int> inv_fact_;

    std::vector<std::vector<int>> block0_;
    std::vector<std::vector<int>> block1_;
    std::vector<Block> triple_blocks_;
    std::vector<Block> pair_blocks_;

    std::vector<std::vector<int>> dp0_;
    std::vector<std::vector<int>> dp1_;

    void build_factorials() {
        fact_.assign(max_blocks_ + 1, 1);
        inv_fact_.assign(max_blocks_ + 1, 1);
        for (int i = 1; i <= max_blocks_; ++i) {
            fact_[i] = static_cast<long long>(fact_[i - 1]) * i % MOD;
        }
        inv_fact_[max_blocks_] = static_cast<int>(mod_pow(fact_[max_blocks_], MOD - 2));
        for (int i = max_blocks_; i >= 1; --i) {
            inv_fact_[i - 1] = static_cast<long long>(inv_fact_[i]) * i % MOD;
        }
    }

    int comb_small(long long N, int k) const {
        if (k < 0 || N < k) return 0;
        long long res = 1;
        for (int i = 0; i < k; ++i) {
            long long term = (N - i) % MOD;
            if (term < 0) term += MOD;
            res = res * term % MOD;
        }
        res = res * inv_fact_[k] % MOD;
        return static_cast<int>(res);
    }

    void build_block_counts() {
        const int S = static_cast<int>(dfa_.trans.size());
        const int size = S * stride_;
        std::vector<int> cur(size, 0);
        std::vector<int> next(size, 0);

        block0_.assign(max_len_ + 1, std::vector<int>(max_tiles_ + 1, 0));
        block1_.assign(max_len_ + 1, std::vector<int>(max_tiles_ + 1, 0));

        cur[0] = 1;

        for (int L = 1; L <= max_len_; ++L) {
            std::fill(next.begin(), next.end(), 0);
            for (int s = 0; s < S; ++s) {
                const int base = s * stride_;
                for (int m = 0; m <= max_tiles_; ++m) {
                    const int val = cur[base + m];
                    if (!val) continue;
                    for (int c = 1; c <= 4; ++c) {
                        const int nm = m + c;
                        if (nm > max_tiles_) break;
                        const int ns = dfa_.trans[s][c];
                        const int idx = ns * stride_ + nm;
                        int v = next[idx] + val;
                        if (v >= MOD) v -= MOD;
                        next[idx] = v;
                    }
                }
            }

            std::vector<int> sum0(max_tiles_ + 1, 0);
            std::vector<int> sum1(max_tiles_ + 1, 0);

            for (int s : dfa_.accept0) {
                const int base = s * stride_;
                for (int m = 0; m <= max_tiles_; ++m) {
                    int v = sum0[m] + next[base + m];
                    if (v >= MOD) v -= MOD;
                    sum0[m] = v;
                }
            }
            for (int s : dfa_.accept1) {
                const int base = s * stride_;
                for (int m = 0; m <= max_tiles_; ++m) {
                    int v = sum1[m] + next[base + m];
                    if (v >= MOD) v -= MOD;
                    sum1[m] = v;
                }
            }
            block0_[L] = std::move(sum0);
            block1_[L] = std::move(sum1);

            cur.swap(next);
        }

        for (int L = 1; L <= max_len_; ++L) {
            for (int M = 0; M <= max_tiles_; ++M) {
                const int c0 = block0_[L][M];
                if (c0 && (M % 3 == 0)) {
                    triple_blocks_.push_back({L, M, c0});
                }
                const int c1 = block1_[L][M];
                if (c1 && (M % 3 == 2)) {
                    pair_blocks_.push_back({L, M, c1});
                }
            }
        }
    }

    void convolve_add(const std::vector<int>& dp,
                      const std::vector<Block>& blocks,
                      std::vector<int>& out) const {
        if (blocks.empty()) return;
        if (threads_ <= 1 || blocks.size() < 200) {
            for (const auto& b : blocks) {
                const int Lb = b.len;
                const int Mb = b.tiles;
                const int cnt = b.count;
                const int Llimit = max_len_ - Lb;
                const int Mlimit = max_tiles_ - Mb;
                for (int L = 0; L <= Llimit; ++L) {
                    const int base = L * stride_;
                    const int out_base = (L + Lb) * stride_ + Mb;
                    for (int M = 0; M <= Mlimit; ++M) {
                        const int val = dp[base + M];
                        if (!val) continue;
                        const int idx = out_base + M;
                        out[idx] = (out[idx] + static_cast<long long>(val) * cnt) % MOD;
                    }
                }
            }
            return;
        }

        const size_t total = out.size();
        const unsigned tcount = std::min<unsigned>(threads_, blocks.size());
        std::vector<std::vector<int>> locals(tcount, std::vector<int>(total, 0));
        std::vector<std::thread> workers;
        workers.reserve(tcount);

        for (unsigned t = 0; t < tcount; ++t) {
            const size_t start = t * blocks.size() / tcount;
            const size_t end = (t + 1) * blocks.size() / tcount;
            workers.emplace_back([&, start, end, t]() {
                auto& local = locals[t];
                for (size_t i = start; i < end; ++i) {
                    const auto& b = blocks[i];
                    const int Lb = b.len;
                    const int Mb = b.tiles;
                    const int cnt = b.count;
                    const int Llimit = max_len_ - Lb;
                    const int Mlimit = max_tiles_ - Mb;
                    for (int L = 0; L <= Llimit; ++L) {
                        const int base = L * stride_;
                        const int out_base = (L + Lb) * stride_ + Mb;
                        for (int M = 0; M <= Mlimit; ++M) {
                            const int val = dp[base + M];
                            if (!val) continue;
                            const int idx = out_base + M;
                            local[idx] = (local[idx] + static_cast<long long>(val) * cnt) % MOD;
                        }
                    }
                }
            });
        }
        for (auto& th : workers) th.join();

        for (size_t i = 0; i < total; ++i) {
            long long sum = out[i];
            for (unsigned t = 0; t < tcount; ++t) {
                sum += locals[t][i];
            }
            out[i] = static_cast<int>(sum % MOD);
        }
    }

    void build_block_dp() {
        const int size = (max_len_ + 1) * stride_;
        dp0_.assign(max_blocks_ + 1, std::vector<int>(size, 0));
        dp1_.assign(max_blocks_ + 1, std::vector<int>(size, 0));
        dp0_[0][0] = 1;

        for (int k = 0; k < max_blocks_; ++k) {
            std::vector<int> next0(size, 0);
            std::vector<int> next1(size, 0);
            convolve_add(dp0_[k], triple_blocks_, next0);
            convolve_add(dp1_[k], triple_blocks_, next1);
            convolve_add(dp0_[k], pair_blocks_, next1);
            dp0_[k + 1] = std::move(next0);
            dp1_[k + 1] = std::move(next1);
        }
    }

    static std::vector<int> poly_mul(const std::vector<int>& a,
                                     const std::vector<int>& b,
                                     int t) {
        std::vector<int> res(t + 1, 0);
        const int n = std::min<int>(t, static_cast<int>(a.size()) - 1);
        const int m = std::min<int>(t, static_cast<int>(b.size()) - 1);
        for (int i = 0; i <= n; ++i) {
            if (!a[i]) continue;
            for (int j = 0; j <= m && i + j <= t; ++j) {
                if (!b[j]) continue;
                res[i + j] = (res[i + j] + static_cast<long long>(a[i]) * b[j]) % MOD;
            }
        }
        return res;
    }

    static std::vector<int> poly_pow(std::vector<int> base, long long exp, int t) {
        std::vector<int> res(t + 1, 0);
        res[0] = 1;
        while (exp > 0) {
            if (exp & 1LL) res = poly_mul(res, base, t);
            exp >>= 1LL;
            if (exp) base = poly_mul(base, base, t);
        }
        return res;
    }
};

}  // namespace

int main() {
    const int max_t = 30;
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    threads = std::min(threads, 8u);

    MahjongCounter counter(max_t, threads);

    bool ok = true;
    auto check = [&](long long n, long long s, int t, long long expected) {
        const long long got = counter.solve(n, s, t);
        if (got != expected) {
            std::cerr << "[CHECK FAILED] w(" << n << ", " << s << ", " << t << ") = "
                      << got << ", expected " << expected << "\n";
            ok = false;
        }
    };

    check(4, 1, 1, 20);
    check(9, 1, 4, 13259);
    check(9, 3, 4, 5237550);
    check(1000, 1000, 5, 107662178);

    if (!ok) return 1;

    const long long n = 100000000LL;
    const long long s = 100000000LL;
    const int t = 30;
    const long long ans = counter.solve(n, s, t);
    std::cout << ans % MOD << "\n";
    return 0;
}
