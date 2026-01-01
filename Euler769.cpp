#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

uint64_t isqrt_u64(uint64_t x) {
    long double r = sqrtl(static_cast<long double>(x));
    uint64_t y = static_cast<uint64_t>(r);
    while ((y + 1) * (y + 1) <= x) ++y;
    while (y * y > x) --y;
    return y;
}

vector<int> build_spf(int limit) {
    vector<int> spf(limit + 1, 0);
    if (limit >= 1) spf[1] = 1;
    for (int i = 2; i <= limit; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            if (1LL * i * i <= limit) {
                for (int j = i * i; j <= limit; j += i) {
                    if (spf[j] == 0) spf[j] = i;
                }
            }
        }
    }
    return spf;
}

struct Counter {
    const vector<int>& spf;
    long double alpha;
    long double beta;
    int inv13[13];

    explicit Counter(const vector<int>& spf_ref) : spf(spf_ref) {
        alpha = 2.0L - sqrtl(3.0L);
        beta = 2.0L + sqrtl(3.0L);
        inv13[0] = 0;
        for (int i = 1; i < 13; ++i) {
            for (int j = 1; j < 13; ++j) {
                if ((i * j) % 13 == 1) {
                    inv13[i] = j;
                    break;
                }
            }
        }
    }

    void build_squarefree_divs(int q, vector<int>& divs, vector<int>& mus) const {
        divs.clear();
        mus.clear();
        divs.push_back(1);
        mus.push_back(1);
        int x = q;
        while (x > 1) {
            int p = spf[x];
            while (x % p == 0) x /= p;
            int current = static_cast<int>(divs.size());
            for (int i = 0; i < current; ++i) {
                divs.push_back(divs[i] * p);
                mus.push_back(-mus[i]);
            }
        }
    }

    long long coprime_count(const vector<int>& divs,
                            const vector<int>& mus,
                            long long k_max) const {
        if (k_max <= 0) return 0;
        long long total = 0;
        for (size_t i = 0; i < divs.size(); ++i) {
            total += static_cast<long long>(mus[i]) * (k_max / divs[i]);
        }
        return total;
    }

    long long coprime_count_residue(const vector<int>& divs,
                                    const vector<int>& mus,
                                    long long k_max,
                                    int residue) const {
        if (k_max <= 0) return 0;
        long long total = 0;
        for (size_t i = 0; i < divs.size(); ++i) {
            int d = divs[i];
            int d_mod = d % 13;
            int inv = inv13[d_mod];
            int t0 = (residue * inv) % 13;
            if (t0 == 0) continue;
            long long first = 1LL * d * t0;
            if (first > k_max) continue;
            total += static_cast<long long>(mus[i]) *
                     (1 + (k_max - first) / (13LL * d));
        }
        return total;
    }

    long long compute(uint64_t N, int threads) const {
        uint64_t q_limit = isqrt_u64(N);
        if (q_limit == 0) return 0;
        if (threads <= 0) {
            unsigned hw = thread::hardware_concurrency();
            threads = hw == 0 ? 1 : static_cast<int>(hw);
        }
        if (q_limit < static_cast<uint64_t>(threads)) {
            threads = static_cast<int>(q_limit);
        }
        vector<long long> partial(threads, 0);
        uint64_t chunk = (q_limit + threads - 1) / threads;

        auto worker = [&](int idx, uint64_t start, uint64_t end) {
            vector<int> divs;
            vector<int> mus;
            divs.reserve(64);
            mus.reserve(64);
            long long local = 0;
            for (uint64_t q = start; q <= end; ++q) {
                build_squarefree_divs(static_cast<int>(q), divs, mus);
                long long q_ll = static_cast<long long>(q);
                uint64_t q2 = static_cast<uint64_t>(q_ll) * q_ll;
                int q_mod13 = static_cast<int>(q % 13);
                int residue = q_mod13 == 0 ? 0 : (13 - (2 * q_mod13) % 13) % 13;

                auto count_with_residue = [&](long long k_max) -> long long {
                    if (k_max <= 0) return 0;
                    long long base = coprime_count(divs, mus, k_max);
                    if (q_mod13 != 0) {
                        base -= coprime_count_residue(divs, mus, k_max, residue);
                    }
                    return base;
                };

                auto count_exclude13 = [&](long long k_max) -> long long {
                    if (k_max <= 0) return 0;
                    long long base = coprime_count(divs, mus, k_max);
                    if (q_mod13 != 0) {
                        base -= coprime_count(divs, mus, k_max / 13);
                    }
                    return base;
                };

                auto zA = [&](long long k) -> long long {
                    __int128 val = static_cast<__int128>(q2) +
                                   static_cast<__int128>(q_ll) * k -
                                   3 * static_cast<__int128>(k) * k;
                    return static_cast<long long>(val);
                };
                auto zC = [&](long long k) -> long long {
                    __int128 val = 13 * static_cast<__int128>(q2) +
                                   13 * static_cast<__int128>(q_ll) * k +
                                   3 * static_cast<__int128>(k) * k;
                    return static_cast<long long>(val);
                };
                auto zD = [&](long long k) -> long long {
                    __int128 val = 3 * static_cast<__int128>(k) * k -
                                   static_cast<__int128>(q_ll) * k -
                                   static_cast<__int128>(q2);
                    return static_cast<long long>(val);
                };

                // Region A: 1 < p/q < 3 - sqrt(3)
                long long k_max = static_cast<long long>(floorl(alpha * q));
                if (k_max > 0) {
                    while (k_max > 0) {
                        __int128 x_val = static_cast<__int128>(q2) -
                                         4 * static_cast<__int128>(q_ll) * k_max +
                                         static_cast<__int128>(k_max) * k_max;
                        if (x_val > 0) break;
                        --k_max;
                    }
                    if (k_max > 0) {
                        long long base = count_with_residue(k_max);
                        __int128 disc128 = 13 * static_cast<__int128>(q2) -
                                           12 * static_cast<__int128>(N);
                        if (disc128 > 0) {
                            uint64_t disc = static_cast<uint64_t>(disc128);
                            uint64_t s = isqrt_u64(disc);
                            long long L = (q_ll - static_cast<long long>(s) + 5) / 6;
                            long long R = (q_ll + static_cast<long long>(s)) / 6;
                            while (L <= k_max && zA(L) <= static_cast<long long>(N)) ++L;
                            while (R >= 1 && zA(R) <= static_cast<long long>(N)) --R;
                            if (L <= R && L <= k_max && R >= 1) {
                                long long L2 = max(1LL, L);
                                long long R2 = min(k_max, R);
                                if (L2 <= R2) {
                                    long long forbidden = count_with_residue(R2) -
                                                          count_with_residue(L2 - 1);
                                    base -= forbidden;
                                }
                            }
                        }
                        local += base;
                    }
                }

                // Common sqrt for regions C and D.
                __int128 disc_plus128 = 13 * static_cast<__int128>(q2) +
                                        12 * static_cast<__int128>(N);
                uint64_t disc_plus = static_cast<uint64_t>(disc_plus128);
                uint64_t s_plus = isqrt_u64(disc_plus);

                // Region C: p/q < -1
                long long k_max_c = (static_cast<long long>(s_plus) - 13LL * q_ll) / 6;
                if (k_max_c > 0) {
                    while (k_max_c > 0 && zC(k_max_c) > static_cast<long long>(N)) {
                        --k_max_c;
                    }
                    while (zC(k_max_c + 1) <= static_cast<long long>(N)) {
                        ++k_max_c;
                    }
                    if (k_max_c > 0) {
                        local += count_exclude13(k_max_c);
                    }
                }

                // Region D: p/q > 3 + sqrt(3)
                long long k_min_d = static_cast<long long>(floorl(beta * q)) + 1;
                if (k_min_d < 1) k_min_d = 1;
                while (true) {
                    __int128 x_val = static_cast<__int128>(q2) -
                                     4 * static_cast<__int128>(q_ll) * k_min_d +
                                     static_cast<__int128>(k_min_d) * k_min_d;
                    if (x_val > 0) break;
                    ++k_min_d;
                }
                long long k_max_d = (static_cast<long long>(s_plus) + q_ll) / 6;
                while (k_max_d > 0 && zD(k_max_d) > static_cast<long long>(N)) {
                    --k_max_d;
                }
                while (zD(k_max_d + 1) <= static_cast<long long>(N)) {
                    ++k_max_d;
                }
                if (k_min_d <= k_max_d) {
                    local += count_with_residue(k_max_d) -
                             count_with_residue(k_min_d - 1);
                }
            }
            partial[idx] = local;
        };

        vector<thread> workers;
        workers.reserve(threads);
        for (int t = 0; t < threads; ++t) {
            uint64_t start = t * chunk + 1;
            uint64_t end = min(q_limit, (t + 1) * chunk);
            if (start > end) {
                partial[t] = 0;
                continue;
            }
            workers.emplace_back(worker, t, start, end);
        }
        for (auto& th : workers) th.join();

        long long total = 0;
        for (long long v : partial) total += v;
        if (N >= 3) total += 1; // (1,1,3)
        return total;
    }
};

} // namespace

int main() {
    const uint64_t N = 100000000000000ULL;
    uint64_t q_limit = isqrt_u64(N);
    vector<int> spf = build_spf(static_cast<int>(q_limit));
    Counter counter(spf);

    // Validation checkpoints.
    long long check1 = counter.compute(1000, 1);
    long long check2 = counter.compute(1000000, 1);
    if (check1 != 142 || check2 != 142463) {
        cerr << "Validation failed: C(1e3)=" << check1
             << " C(1e6)=" << check2 << "\n";
        return 1;
    }

    int threads = static_cast<int>(thread::hardware_concurrency());
    long long answer = counter.compute(N, threads);
    cout << answer << "\n";
    return 0;
}
