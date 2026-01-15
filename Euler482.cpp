#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using int64 = long long;
using u128 = unsigned __int128;

struct Leg {
    int64 x;
    int64 t;
};

int max_inradius_for_perimeter(int64 P) {
    long double bound = static_cast<long double>(P) * std::sqrt(3.0L) / 18.0L;
    int r_max = static_cast<int>(std::floor(bound + 1e-12L));
    return std::max(0, r_max);
}

std::vector<int> build_spf(int limit) {
    std::vector<int> spf(limit + 1, 0);
    if (limit >= 1) spf[1] = 1;
    for (int i = 2; i <= limit; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            if ((int64)i * i <= limit) {
                for (int64 j = (int64)i * i; j <= limit; j += i) {
                    if (spf[(int)j] == 0) spf[(int)j] = i;
                }
            }
        }
    }
    return spf;
}

void factorize(int n, const std::vector<int>& spf, std::vector<std::pair<int, int>>& out) {
    out.clear();
    while (n > 1) {
        int p = spf[n];
        int cnt = 0;
        while (n % p == 0) {
            n /= p;
            ++cnt;
        }
        out.push_back({p, cnt});
    }
}

std::vector<Leg> build_legs(int r, int64 r2, int64 max_x, const std::vector<int>& spf) {
    std::vector<std::pair<int, int>> factors;
    factorize(r, spf, factors);

    std::vector<int64> divisors;
    divisors.reserve(64);
    divisors.push_back(1);
    for (const auto& f : factors) {
        int p = f.first;
        int exp = 2 * f.second;
        int64 p_pow = 1;
        const std::size_t base = divisors.size();
        for (int e = 1; e <= exp; ++e) {
            p_pow *= p;
            for (std::size_t i = 0; i < base; ++i) {
                divisors.push_back(divisors[i] * p_pow);
            }
        }
    }

    std::vector<Leg> legs;
    legs.reserve(divisors.size());
    for (int64 d : divisors) {
        if (d > r) continue;
        int64 d2 = r2 / d;
        if (((d + d2) & 1LL) != 0) continue;
        int64 x = (d2 - d) / 2;
        if (x <= 0 || x > max_x) continue;
        int64 t = (d2 + d) / 2;
        legs.push_back({x, t});
    }
    std::sort(legs.begin(), legs.end(), [](const Leg& a, const Leg& b) { return a.x < b.x; });
    return legs;
}

u128 compute_sum_range(int64 P, int r_start, int r_end, const std::vector<int>& spf) {
    const int64 P_half = P / 2;
    u128 total = 0;
    for (int r = r_start; r <= r_end; ++r) {
        int64 r2 = (int64)r * r;
        auto legs = build_legs(r, r2, P_half, spf);
        if (legs.size() < 2) continue;

        std::vector<int64> xs;
        std::vector<int64> ts;
        xs.reserve(legs.size());
        ts.reserve(legs.size());
        for (const auto& leg : legs) {
            xs.push_back(leg.x);
            ts.push_back(leg.t);
        }

        const std::size_t n = xs.size();
        for (std::size_t i = 0; i < n; ++i) {
            int64 x = xs[i];
            int64 tx = ts[i];
            for (std::size_t j = i; j < n; ++j) {
                int64 y = xs[j];
                int64 ty = ts[j];
                int64 denom = x * y - r2;
                if (denom <= 0) continue;
                u128 numer = (u128)r2 * (x + y);
                if (numer % denom != 0) continue;
                u128 z_u = numer / denom;
                int64 max_z = P_half - x - y;
                if (max_z <= 0) continue;
                if (z_u > (u128)max_z) continue;
                int64 z = (int64)z_u;
                if (z < y) continue;
                auto it = std::lower_bound(xs.begin(), xs.end(), z);
                if (it == xs.end() || *it != z) continue;
                std::size_t k = static_cast<std::size_t>(it - xs.begin());
                int64 tz = ts[k];

                int64 p = 2 * (x + y + z);
                int64 L = p + tx + ty + tz;
                total += (u128)L;
            }
        }
    }
    return total;
}

u128 compute_sum(int64 P, const std::vector<int>& spf, int r_max, unsigned threads) {
    if (r_max <= 0) return 0;
    if (threads == 0) threads = 1;
    if (threads == 1 || r_max < 200) {
        return compute_sum_range(P, 1, r_max, spf);
    }

    int chunk = (r_max + (int)threads - 1) / (int)threads;
    std::vector<std::thread> pool;
    std::vector<u128> partial(threads, 0);
    pool.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
        int start = (int)t * chunk + 1;
        int end = std::min(r_max, start + chunk - 1);
        if (start > end) continue;
        pool.emplace_back([&, start, end, t]() {
            partial[t] = compute_sum_range(P, start, end, spf);
        });
    }
    for (auto& th : pool) th.join();

    u128 total = 0;
    for (const auto& v : partial) total += v;
    return total;
}

std::string to_string_u128(u128 value) {
    if (value == 0) return "0";
    std::string out;
    while (value > 0) {
        int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

int64 brute_sum(int P) {
    int64 total = 0;
    for (int a = 1; a <= P / 3; ++a) {
        for (int b = a; b <= (P - a) / 2; ++b) {
            int max_c = std::min(P - a - b, a + b - 1);
            if (max_c < b) continue;
            for (int c = b; c <= max_c; ++c) {
                int p = a + b + c;
                if (p > P) break;
                double s = 0.5 * p;
                double x = s - a;
                double y = s - b;
                double z = s - c;
                double area2 = s * x * y * z;
                if (area2 <= 0.0) continue;
                double r = std::sqrt(area2) / s;
                double IA = std::sqrt(r * r + x * x);
                double IB = std::sqrt(r * r + y * y);
                double IC = std::sqrt(r * r + z * z);
                auto round_int = [](double v) { return std::llround(v); };
                long long ia = round_int(IA);
                long long ib = round_int(IB);
                long long ic = round_int(IC);
                if (std::fabs(IA - ia) > 1e-9) continue;
                if (std::fabs(IB - ib) > 1e-9) continue;
                if (std::fabs(IC - ic) > 1e-9) continue;
                total += p + ia + ib + ic;
            }
        }
    }
    return total;
}

} // namespace

int main(int argc, char** argv) {
    int64 P = 10000000LL;
    if (argc >= 2) {
        P = std::stoll(argv[1]);
    }
    unsigned threads = std::thread::hardware_concurrency();
    if (argc >= 3) {
        threads = static_cast<unsigned>(std::stoul(argv[2]));
    }

    int r_max = max_inradius_for_perimeter(P);
    auto spf = build_spf(r_max);

    if (P <= 200) {
        u128 fast = compute_sum(P, spf, r_max, 1);
        int64 slow = brute_sum(static_cast<int>(P));
        if (fast != static_cast<u128>(slow)) {
            std::cerr << "Validation failed for P=" << P << ": fast="
                      << to_string_u128(fast) << ", brute=" << slow << "\n";
            return 1;
        }
        std::cout << to_string_u128(fast) << "\n";
        return 0;
    }

    if (P >= 1000) {
        int r_max_check = max_inradius_for_perimeter(1000);
        u128 check = compute_sum(1000, spf, r_max_check, 1);
        if (check != 3619) {
            std::cerr << "Validation failed for P=1000: got "
                      << to_string_u128(check) << ", expected 3619\n";
            return 1;
        }
    }

    u128 answer = compute_sum(P, spf, r_max, threads);
    std::cout << to_string_u128(answer) << "\n";
    return 0;
}
