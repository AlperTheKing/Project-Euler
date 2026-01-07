#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

namespace {

constexpr int64_t MOD = 912491249LL;

inline int64_t mod_add(int64_t a, int64_t b) {
    int64_t v = a + b;
    if (v >= MOD) v -= MOD;
    return v;
}

inline int64_t mod_mul(int64_t a, int64_t b) {
    return static_cast<int64_t>((__int128)a * b % MOD);
}

struct PeriodInfo {
    int period = 0;
    int start = 0;
    vector<int> seq;
};

struct GrundyData {
    vector<int> g;
    int max_g = 0;
    PeriodInfo even;
    PeriodInfo odd;
};

vector<int> compute_grundy(int max_n, int& max_g) {
    vector<int> g(static_cast<size_t>(max_n + 1), 0);
    max_g = 0;
    const int moves[] = {1, 2, 4, 9};

    for (int n = 1; n <= max_n; ++n) {
        uint64_t mask = 0;
        for (int mv : moves) {
            if (n >= mv) {
                int val = g[n - mv];
                if (val >= 63) {
                    cerr << "Grundy value too large for mask\n";
                    exit(1);
                }
                mask |= (1ULL << val);
            }
        }
        int half = n / 2;
        for (int a = 1; a <= half; ++a) {
            int val = g[a] ^ g[n - a];
            if (val >= 63) {
                cerr << "Grundy value too large for mask\n";
                exit(1);
            }
            mask |= (1ULL << val);
        }
        int mex = 0;
        while (mask & (1ULL << mex)) ++mex;
        if (mex >= 63) {
            cerr << "Grundy value too large for mask\n";
            exit(1);
        }
        g[n] = mex;
        if (mex > max_g) max_g = mex;
    }
    return g;
}

PeriodInfo detect_period(const vector<int>& seq, int max_period = 300, int min_reps = 6) {
    const int n = static_cast<int>(seq.size());
    for (int p = 1; p <= max_period; ++p) {
        int L = p * min_reps;
        if (L + p > n) continue;
        int start_tail = n - L;
        bool ok = true;
        for (int i = start_tail; i + p < n; ++i) {
            if (seq[i] != seq[i + p]) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;
        for (int s = 0; s + p < n; ++s) {
            bool good = true;
            for (int i = s; i + p < n; ++i) {
                if (seq[i] != seq[i + p]) {
                    good = false;
                    break;
                }
            }
            if (good) {
                PeriodInfo info;
                info.period = p;
                info.start = s;
                info.seq = seq;
                return info;
            }
        }
    }
    return PeriodInfo{};
}

GrundyData build_grundy_data(int max_n) {
    GrundyData data;
    data.g = compute_grundy(max_n, data.max_g);

    vector<int> even;
    vector<int> odd;
    even.reserve(static_cast<size_t>(max_n / 2));
    odd.reserve(static_cast<size_t>((max_n + 1) / 2));

    for (int k = 1; 2 * k <= max_n; ++k) {
        even.push_back(data.g[2 * k]);
    }
    for (int k = 1; 2 * k - 1 <= max_n; ++k) {
        odd.push_back(data.g[2 * k - 1]);
    }

    data.even = detect_period(even);
    data.odd = detect_period(odd);
    if (data.even.period == 0 || data.odd.period == 0) {
        cerr << "Failed to detect periods\n";
        exit(1);
    }
    return data;
}

void add_counts(const PeriodInfo& info, long long total, vector<long long>& counts) {
    if (total <= 0) return;
    const vector<int>& seq = info.seq;
    const int start = info.start;
    const int period = info.period;

    if (total <= start) {
        for (long long i = 0; i < total; ++i) {
            counts[seq[static_cast<size_t>(i)]]++;
        }
        return;
    }

    for (int i = 0; i < start; ++i) {
        counts[seq[static_cast<size_t>(i)]]++;
    }

    long long rem = total - start;
    vector<long long> freq(counts.size(), 0);
    for (int i = 0; i < period; ++i) {
        freq[seq[static_cast<size_t>(start + i)]]++;
    }
    long long cycles = rem / period;
    long long tail = rem % period;
    for (size_t g = 0; g < counts.size(); ++g) {
        counts[g] += freq[g] * cycles;
    }
    for (int i = 0; i < tail; ++i) {
        counts[seq[static_cast<size_t>(start + i)]]++;
    }
}

vector<long long> grundy_counts(long long N, const GrundyData& data) {
    vector<long long> counts(static_cast<size_t>(data.max_g + 1), 0);
    long long total_even = N / 2;
    long long total_odd = (N + 1) / 2;
    add_counts(data.even, total_even, counts);
    add_counts(data.odd, total_odd, counts);
    return counts;
}

vector<int64_t> build_comb(long long c, int m) {
    vector<int64_t> comb(static_cast<size_t>(m + 1), 0);
    cpp_int val = 1;
    comb[0] = 1;
    for (int k = 1; k <= m; ++k) {
        val *= static_cast<long long>(c + k - 1);
        val /= k;
        comb[static_cast<size_t>(k)] = (val % MOD).convert_to<long long>();
    }
    return comb;
}

vector<vector<int64_t>> update_dp(const vector<vector<int64_t>>& dp,
                                  int gval,
                                  const vector<pair<int, int64_t>>& even_terms,
                                  const vector<pair<int, int64_t>>& odd_terms,
                                  int m) {
    const int states = static_cast<int>(dp.size());
    int thread_count = static_cast<int>(thread::hardware_concurrency());
    if (thread_count <= 0) thread_count = 1;
    if (thread_count > states) thread_count = states;

    vector<vector<int64_t>> next(states, vector<int64_t>(static_cast<size_t>(m + 1), 0));
    if (thread_count == 1) {
        for (int mask = 0; mask < states; ++mask) {
            const auto& cur = dp[mask];
            auto& dest_same = next[mask];
            auto& dest_xor = next[mask ^ gval];
            for (int i = 0; i <= m; ++i) {
                int64_t curv = cur[static_cast<size_t>(i)];
                if (curv == 0) continue;
                for (const auto& term : even_terms) {
                    int j = term.first;
                    if (i + j > m) break;
                    int64_t add = mod_mul(curv, term.second);
                    dest_same[static_cast<size_t>(i + j)] =
                        mod_add(dest_same[static_cast<size_t>(i + j)], add);
                }
                for (const auto& term : odd_terms) {
                    int j = term.first;
                    if (i + j > m) break;
                    int64_t add = mod_mul(curv, term.second);
                    dest_xor[static_cast<size_t>(i + j)] =
                        mod_add(dest_xor[static_cast<size_t>(i + j)], add);
                }
            }
        }
        return next;
    }

    vector<vector<vector<int64_t>>> partial(
        static_cast<size_t>(thread_count),
        vector<vector<int64_t>>(states, vector<int64_t>(static_cast<size_t>(m + 1), 0)));

    vector<thread> workers;
    workers.reserve(static_cast<size_t>(thread_count));
    for (int t = 0; t < thread_count; ++t) {
        workers.emplace_back([&, t]() {
            for (int mask = t; mask < states; mask += thread_count) {
                const auto& cur = dp[mask];
                auto& dest_same = partial[static_cast<size_t>(t)][mask];
                auto& dest_xor = partial[static_cast<size_t>(t)][mask ^ gval];
                for (int i = 0; i <= m; ++i) {
                    int64_t curv = cur[static_cast<size_t>(i)];
                    if (curv == 0) continue;
                    for (const auto& term : even_terms) {
                        int j = term.first;
                        if (i + j > m) break;
                        int64_t add = mod_mul(curv, term.second);
                        dest_same[static_cast<size_t>(i + j)] =
                            mod_add(dest_same[static_cast<size_t>(i + j)], add);
                    }
                    for (const auto& term : odd_terms) {
                        int j = term.first;
                        if (i + j > m) break;
                        int64_t add = mod_mul(curv, term.second);
                        dest_xor[static_cast<size_t>(i + j)] =
                            mod_add(dest_xor[static_cast<size_t>(i + j)], add);
                    }
                }
            }
        });
    }
    for (auto& th : workers) th.join();

    for (int t = 0; t < thread_count; ++t) {
        for (int mask = 0; mask < states; ++mask) {
            auto& dest = next[mask];
            const auto& src = partial[static_cast<size_t>(t)][mask];
            for (int i = 0; i <= m; ++i) {
                int64_t val = src[static_cast<size_t>(i)];
                if (val == 0) continue;
                dest[static_cast<size_t>(i)] = mod_add(dest[static_cast<size_t>(i)], val);
            }
        }
    }

    return next;
}

int64_t compute_S(long long N, int m, const GrundyData& data) {
    vector<long long> counts = grundy_counts(N, data);

    int max_g = data.max_g;
    int bits = 0;
    while ((1 << bits) <= max_g) ++bits;
    int states = 1 << bits;

    vector<vector<int64_t>> dp(states, vector<int64_t>(static_cast<size_t>(m + 1), 0));
    dp[0][0] = 1;

    for (int gval = 0; gval <= max_g; ++gval) {
        long long c = counts[static_cast<size_t>(gval)];
        if (c == 0) continue;

        vector<int64_t> comb = build_comb(c, m);
        vector<pair<int, int64_t>> even_terms;
        vector<pair<int, int64_t>> odd_terms;
        even_terms.reserve(static_cast<size_t>(m / 2 + 1));
        odd_terms.reserve(static_cast<size_t>(m / 2 + 1));
        for (int k = 0; k <= m; ++k) {
            int64_t v = comb[static_cast<size_t>(k)];
            if (v == 0) continue;
            if (k & 1) {
                odd_terms.push_back({k, v});
            } else {
                even_terms.push_back({k, v});
            }
        }

        dp = update_dp(dp, gval, even_terms, odd_terms, m);
    }

    return dp[0][static_cast<size_t>(m)];
}

}  // namespace

int main() {
    const int max_n = 20000;
    GrundyData data = build_grundy_data(max_n);

    const int64_t check1 = compute_S(12, 4, data);
    if (check1 != 204) {
        cerr << "Validation failure: S(12,4)\n";
        return 1;
    }

    const int64_t check2 = compute_S(124, 9, data);
    const int64_t expected2 = 2259208528408LL % MOD;
    if (check2 != expected2) {
        cerr << "Validation failure: S(124,9)\n";
        return 1;
    }

    const long long N = 12'491'249LL;
    const int m = 1249;
    cout << compute_S(N, m, data) << '\n';
    return 0;
}
