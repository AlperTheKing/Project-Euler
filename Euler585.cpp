#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

using int64 = long long;

static std::vector<int> totient_sieve(int n) {
    std::vector<int> phi(n + 1);
    std::iota(phi.begin(), phi.end(), 0);
    for (int i = 2; i <= n; ++i) {
        if (phi[i] == i) {
            for (int j = i; j <= n; j += i) {
                phi[j] -= phi[j] / i;
            }
        }
    }
    return phi;
}

// f(s) = phi(s)/2 - R(s), where R(s) counts primitive i^2 + j^2 representations.
static std::vector<int> build_f(int n) {
    std::vector<int> phi = totient_sieve(n);
    std::vector<int> f(n + 1, 0);
    for (int s = 3; s <= n; ++s) {
        f[s] = phi[s] / 2;
    }

    int lim = static_cast<int>(std::sqrt(static_cast<long double>(n)));
    for (int i = 2; i <= lim; ++i) {
        int i2 = i * i;
        for (int j = 1; j < i; ++j) {
            if (std::gcd(i, j) != 1) continue;
            int s = i2 + j * j;
            if (s > n) break;
            f[s] -= 1;
        }
    }
    return f;
}

static std::vector<int64> prefix_sum(const std::vector<int>& f) {
    std::vector<int64> pref(f.size(), 0);
    for (size_t i = 1; i < f.size(); ++i) {
        pref[i] = pref[i - 1] + f[i];
    }
    return pref;
}

static int64 partial_sum_floor(int64 m, const std::vector<int64>& pref) {
    int64 ret = 0;
    for (int64 i = 1; i <= m;) {
        int64 q = m / i;
        int64 j = m / q + 1;
        ret += q * (pref[j - 1] - pref[i - 1]);
        i = j;
    }
    return ret;
}

static int64 count_type1(int n, const std::vector<int>& f) {
    int64 cnt = 0;
    for (int s = 3; s <= n; ++s) {
        if (f[s] == 0) continue;
        cnt += static_cast<int64>(f[s]) * (n / s);
    }
    return cnt;
}

static int64 count_type2_ordered(int n, const std::vector<int64>& pref) {
    int64 cnt = 0;
    for (int64 i = 1; i <= n;) {
        int64 q = n / i;
        int64 j = n / q + 1;
        int64 sum_f = pref[j - 1] - pref[i - 1];
        cnt += sum_f * partial_sum_floor(q, pref);
        i = j;
    }
    return cnt;
}

static std::vector<char> squarefree_flags(int n) {
    std::vector<char> ok(n + 1, 1);
    ok[0] = 0;
    for (int p = 2; static_cast<int64>(p) * p <= n; ++p) {
        int64 p2 = static_cast<int64>(p) * p;
        for (int64 j = p2; j <= n; j += p2) {
            ok[static_cast<size_t>(j)] = 0;
        }
    }
    return ok;
}

static int64 count_overlap(int n) {
    int lim = static_cast<int>(std::sqrt(static_cast<long double>(n)));
    auto is_sf = squarefree_flags(lim);

    std::vector<int> sf;
    sf.reserve(lim);
    for (int i = 1; i <= lim; ++i) {
        if (is_sf[i]) sf.push_back(i);
    }

    std::vector<int64> squares(lim + 1, 0);
    for (int i = 0; i <= lim; ++i) squares[i] = static_cast<int64>(i) * i;

    auto worker = [&](int start, int end) {
        int64 acc = 0;
        for (int idx = start; idx < end; ++idx) {
            int p = sf[idx];
            if (static_cast<int64>(p + 1) * (p + 1) > n) break;
            for (int q : sf) {
                if (std::gcd(p, q) != 1) continue;
                int64 pq = static_cast<int64>(p) * q;
                if ((pq + 1) * (p + q) > n) break;
                for (int r : sf) {
                    if (std::gcd(pq, static_cast<int64>(r)) != 1) continue;
                    int64 pr = static_cast<int64>(p) * r;
                    if ((pq + r) * (pr + q) > n) break;
                    int64 pqr = pq * r;
                    for (int s : sf) {
                        if (std::gcd(pqr, static_cast<int64>(s)) != 1) continue;
                        int64 rs = static_cast<int64>(r) * s;
                        int64 qs = static_cast<int64>(q) * s;
                        if ((pq + rs) * (pr + qs) > n) break;

                        int64 u = pq;
                        int64 v = rs;
                        int64 a = pr;
                        int64 b = qs;
                        if (u == v || a == b) continue;

                        int64 max_x = n / (a + b);
                        for (int w1 = 1; w1 <= lim; ++w1) {
                            int64 uw1 = u * squares[w1];
                            if (uw1 + v > max_x) break;
                            for (int w2 = 1; w2 <= lim; ++w2) {
                                int64 vw2 = v * squares[w2];
                                int64 x = uw1 + vw2;
                                if (x > max_x) break;
                                if (uw1 <= vw2) continue;
                                if (std::gcd(uw1, vw2) != 1) continue;

                                int64 m = n / x;
                                for (int w3 = 1; w3 <= lim; ++w3) {
                                    int64 aw3 = a * squares[w3];
                                    if (aw3 + b > m) break;
                                    for (int w4 = 1; w4 <= lim; ++w4) {
                                        int64 bw4 = b * squares[w4];
                                        int64 y = aw3 + bw4;
                                        if (y > m) break;
                                        if (aw3 <= bw4) continue;
                                        if (std::gcd(aw3, bw4) != 1) continue;
                                        acc += m / y;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return acc;
    };

    unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    if (threads == 1 || n < 200000) {
        return worker(0, static_cast<int>(sf.size()));
    }

    std::vector<int64> partial(threads, 0);
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) {
        int start = static_cast<int>((static_cast<int64>(sf.size()) * t) / threads);
        int end = static_cast<int>((static_cast<int64>(sf.size()) * (t + 1)) / threads);
        pool.emplace_back([&, t, start, end]() {
            partial[t] = worker(start, end);
        });
    }
    for (auto& th : pool) th.join();

    int64 total = 0;
    for (int64 v : partial) total += v;
    return total;
}

static int64 compute_f_value(int n) {
    auto f = build_f(n);
    auto pref = prefix_sum(f);

    int64 cnt2 = count_type1(n, f);
    int64 cnt1 = count_type2_ordered(n, pref);
    int64 cnt3 = count_overlap(n);

    return cnt2 + (cnt1 - cnt3) / 2;
}

int main(int argc, char** argv) {
    int n = 5000000;
    if (argc > 1) n = std::atoi(argv[1]);

    std::vector<std::pair<int, int64>> checks = {
        {10, 17}, {15, 46}, {20, 86}, {30, 213}, {100, 2918}, {5000, 11134074}
    };
    for (const auto& c : checks) {
        int64 got = compute_f_value(c.first);
        if (got != c.second) {
            std::cerr << "Validation failed at n=" << c.first << ": got=" << got
                      << " expected=" << c.second << "\n";
            return 1;
        }
    }

    std::cout << compute_f_value(n) << "\n";
    return 0;
}
