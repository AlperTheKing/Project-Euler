#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

struct Mask {
    uint64_t lo;  // bits 0..63
    uint64_t hi;  // bits 64..80 stored in bits 0..16
};

static inline int popcnt_u64(uint64_t x) {
    return __builtin_popcountll(x);
}

static inline int mask_size(const Mask& m) {
    return popcnt_u64(m.lo) + popcnt_u64(m.hi);
}

static cpp_int binom_int(int n, int k) {
    if (k < 0 || k > n) return 0;
    k = min(k, n - k);
    cpp_int res = 1;
    for (int i = 1; i <= k; ++i) {
        res *= (n - k + i);
        res /= i;
    }
    return res;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int threads = static_cast<int>(thread::hardware_concurrency());
    if (threads <= 0) threads = 4;
    bool validate = true;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--no-validate") {
            validate = false;
        } else if (arg == "--validate") {
            validate = true;
        } else {
            threads = max(1, atoi(arg.c_str()));
        }
    }

    // Points are GF(3)^4, indexed by base-3 expansion.
    array<array<int, 4>, 81> coord{};
    for (int id = 0; id < 81; ++id) {
        int t = id;
        for (int d = 0; d < 4; ++d) {
            coord[id][d] = t % 3;
            t /= 3;
        }
    }
    const int pow3[4] = {1, 3, 9, 27};
    auto idx_from_coord = [&](const array<int, 4>& c) -> int {
        int id = 0;
        for (int d = 0; d < 4; ++d) id += c[d] * pow3[d];
        return id;
    };

    // Generate all affine lines (SETs) via point pairs.
    unordered_set<uint32_t> seen;
    seen.reserve(2000);
    vector<array<int, 3>> lines;
    lines.reserve(1080);

    for (int i = 0; i < 81; ++i) {
        for (int j = i + 1; j < 81; ++j) {
            array<int, 4> c{};
            for (int d = 0; d < 4; ++d) {
                int s = coord[i][d] + coord[j][d];
                c[d] = (3 - (s % 3)) % 3;
            }
            int k = idx_from_coord(c);

            array<int, 3> tri = {i, j, k};
            sort(tri.begin(), tri.end());
            uint32_t key = static_cast<uint32_t>(tri[0] | (tri[1] << 7) | (tri[2] << 14));
            if (seen.insert(key).second) {
                lines.push_back(tri);
            }
        }
    }

    if (validate && static_cast<int>(lines.size()) != 1080) {
        cerr << "ERROR: expected 1080 lines, got " << lines.size() << "\n";
        return 1;
    }

    vector<Mask> masks;
    masks.reserve(lines.size());
    for (const auto& tri : lines) {
        Mask m{0, 0};
        for (int p : tri) {
            if (p < 64) m.lo |= (1ULL << p);
            else m.hi |= (1ULL << (p - 64));
        }
        masks.push_back(m);

        if (validate && mask_size(m) != 3) {
            cerr << "ERROR: line mask size != 3\n";
            return 1;
        }
    }

    if (validate) {
        array<int, 81> point_deg{};
        for (const auto& tri : lines) {
            for (int p : tri) point_deg[p]++;
        }
        for (int p = 0; p < 81; ++p) {
            if (point_deg[p] != 40) {
                cerr << "ERROR: point " << p << " has degree " << point_deg[p]
                     << " (expected 40)\n";
                return 1;
            }
        }
    }

    const int L = static_cast<int>(masks.size());
    const Mask fixed = masks[0];

    vector<array<uint64_t, 13>> thread_counts(threads);
    auto worker = [&](int tid, int l2_begin, int l2_end) {
        array<uint64_t, 13> local{};
        for (int l2 = l2_begin; l2 < l2_end; ++l2) {
            const Mask& m2 = masks[l2];
            const uint64_t u12_lo = fixed.lo | m2.lo;
            const uint64_t u12_hi = fixed.hi | m2.hi;

            for (int l3 = 0; l3 < L; ++l3) {
                const Mask& m3 = masks[l3];
                const uint64_t u123_lo = u12_lo | m3.lo;
                const uint64_t u123_hi = u12_hi | m3.hi;

                for (int l4 = 0; l4 < L; ++l4) {
                    const Mask& m4 = masks[l4];
                    uint64_t u_lo = u123_lo | m4.lo;
                    uint64_t u_hi = u123_hi | m4.hi;
                    int sz = popcnt_u64(u_lo) + popcnt_u64(u_hi);
                    local[sz]++;
                }
            }
        }
        thread_counts[tid] = local;
    };

    vector<thread> pool;
    pool.reserve(threads);
    int chunk = (L + threads - 1) / threads;
    for (int t = 0; t < threads; ++t) {
        int b = t * chunk;
        int e = min(L, b + chunk);
        if (b >= e) {
            thread_counts[t].fill(0);
            continue;
        }
        pool.emplace_back(worker, t, b, e);
    }
    for (auto& th : pool) th.join();

    array<uint64_t, 13> counts_fixed{};
    for (int t = 0; t < threads; ++t) {
        for (int sz = 0; sz <= 12; ++sz) counts_fixed[sz] += thread_counts[t][sz];
    }

    if (validate) {
        const unsigned long long expected_triples = 1ULL * L * L * L;
        unsigned long long got_triples = 0;
        for (int sz = 0; sz <= 12; ++sz) got_triples += counts_fixed[sz];
        if (got_triples != expected_triples) {
            cerr << "ERROR: triple count mismatch: got " << got_triples
                 << ", expected " << expected_triples << "\n";
            return 1;
        }
    }

    array<uint64_t, 13> N{};
    for (int sz = 0; sz <= 12; ++sz) {
        N[sz] = counts_fixed[sz] * static_cast<uint64_t>(L);
    }

    if (validate) {
        const unsigned long long expected_quads = 1ULL * L * L * L * L;
        unsigned long long got_quads = 0;
        for (int sz = 0; sz <= 12; ++sz) got_quads += N[sz];
        if (got_quads != expected_quads) {
            cerr << "ERROR: quad count mismatch\n";
            return 1;
        }
    }

    if (validate) {
        cpp_int lhs = 0;
        for (int sz = 0; sz <= 12; ++sz) lhs += cpp_int(sz) * cpp_int(N[sz]);

        const long long diff = 74465;     // 27^4 - 26^4
        const long long forty4 = 2560000; // 40^4
        cpp_int rhs = cpp_int(81) * cpp_int(diff) * cpp_int(forty4);
        if (lhs != rhs) {
            cerr << "ERROR: expected union-size total mismatch\n";
            return 1;
        }
    }

    auto compute_F = [&](int n) -> cpp_int {
        cpp_int sum = 0;
        for (int m = 0; m <= 12; ++m) {
            if (N[m] == 0 || n < m) continue;
            sum += cpp_int(N[m]) * binom_int(81 - m, n - m);
        }
        return sum;
    };

    if (validate) {
        cpp_int F3 = compute_F(3);
        cpp_int F6 = compute_F(6);
        if (F3 != 1080) {
            cerr << "ERROR: F(3) mismatch: got " << F3 << ", expected 1080\n";
            return 1;
        }
        if (F6 != 159690960) {
            cerr << "ERROR: F(6) mismatch: got " << F6 << ", expected 159690960\n";
            return 1;
        }
    }

    cpp_int F12 = compute_F(12);
    cout << F12 << "\n";
    return 0;
}
