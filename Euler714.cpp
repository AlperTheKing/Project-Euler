#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <limits>
#include <numeric>
#include <pthread.h>
#include <string>
#include <iostream>
#include <unistd.h>
#include <vector>

namespace {

using u128 = unsigned __int128;

u128 gcd_u128(u128 a, u128 b) {
    while (b != 0) {
        const u128 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

std::int64_t mod_inverse(std::int64_t a, std::int64_t mod) {
    std::int64_t b = mod;
    std::int64_t x0 = 1;
    std::int64_t x1 = 0;
    while (b != 0) {
        const std::int64_t q = a / b;
        const std::int64_t t = a % b;
        a = b;
        b = t;
        const std::int64_t nx = x0 - q * x1;
        x0 = x1;
        x1 = nx;
    }
    if (x0 < 0) {
        x0 %= mod;
        x0 += mod;
    }
    return x0 % mod;
}

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        s.push_back(static_cast<char>('0' + static_cast<int>(x % 10)));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::string scientific_13sig_u128(u128 value) {
    std::string digits = to_string_u128(value);
    int exponent = static_cast<int>(digits.size()) - 1;

    std::string sig;
    if (digits.size() >= 13) {
        sig = digits.substr(0, 13);
        const int round_digit = (digits.size() > 13) ? (digits[13] - '0') : 0;
        if (round_digit >= 5) {
            int i = 12;
            while (i >= 0 && sig[static_cast<std::size_t>(i)] == '9') {
                sig[static_cast<std::size_t>(i)] = '0';
                --i;
            }
            if (i >= 0) {
                ++sig[static_cast<std::size_t>(i)];
            } else {
                sig = "1" + sig;
                ++exponent;
            }
        }
        if (sig.size() > 13) {
            sig = sig.substr(0, 13);
        }
    } else {
        sig = digits;
        sig.append(13 - sig.size(), '0');
    }

    return std::string(1, sig[0]) + "." + sig.substr(1) + "e" + std::to_string(exponent);
}

class Solver {
  private:
    struct RangeValue {
        u128 max_sum = 0;
        u128 min_sum = 0;
    };

    struct ResidueSet {
        explicit ResidueSet(int mod) : index(static_cast<std::size_t>(mod), -1), n(mod) {}

        void insert_or_update(int residue, u128 max_sum, u128 min_sum) {
            int& slot = index[static_cast<std::size_t>(residue)];
            if (slot < 0) {
                slot = static_cast<int>(keys.size());
                keys.push_back(residue);
                vals.push_back(RangeValue{max_sum, min_sum});
                return;
            }
            RangeValue& cur = vals[static_cast<std::size_t>(slot)];
            if (max_sum > cur.max_sum) cur.max_sum = max_sum;
            if (min_sum < cur.min_sum) cur.min_sum = min_sum;
        }

        void add_shift(int shift_residue, u128 weight) {
            const std::vector<int> old_keys = keys;
            const std::vector<RangeValue> old_vals = vals;
            for (std::size_t i = 0; i < old_keys.size(); ++i) {
                int nr = old_keys[i] + shift_residue;
                if (nr >= n) nr -= n;
                insert_or_update(nr, old_vals[i].max_sum + weight, old_vals[i].min_sum + weight);
            }
        }

        int lookup_index(int residue) const {
            return index[static_cast<std::size_t>(residue)];
        }

        std::size_t size() const { return keys.size(); }

        std::vector<int> index;
        std::vector<int> keys;
        std::vector<RangeValue> vals;
        int n;
    };

    struct GInfo {
        int g = 0;
        int gcd_ng = 1;
        int ng = 1;
        int inv = 0;
    };

    struct WorkerArgs {
        const Solver* solver = nullptr;
        int k = 0;
        std::atomic<int>* next = nullptr;
        int chunk = 0;
        u128* out = nullptr;
    };

    static void* worker_entry(void* raw) {
        auto* args = static_cast<WorkerArgs*>(raw);
        u128 local = 0;
        while (true) {
            const int start = args->next->fetch_add(args->chunk, std::memory_order_relaxed);
            if (start > args->k) break;
            const int end = std::min(args->k, start + args->chunk - 1);
            for (int n = start; n <= end; ++n) {
                local += args->solver->d_u128(n);
            }
        }
        *args->out = local;
        return nullptr;
    }

  public:
    u128 d_u128(int n, int start = 0) const {
        std::vector<GInfo> g_infos;
        g_infos.reserve(18);
        for (int g = -9; g <= 9; ++g) {
            if (g == 0) continue;
            const int g_abs = (g < 0) ? -g : g;
            const int d = std::gcd(g_abs, n);
            const int ng = n / d;
            int g_reduced = g / d;
            g_reduced %= ng;
            if (g_reduced < 0) g_reduced += ng;
            const int inv = static_cast<int>(mod_inverse(g_reduced, ng));
            g_infos.push_back(GInfo{g, d, ng, inv});
        }

        ResidueSet xset(n);
        ResidueSet yset(n);
        xset.insert_or_update(0, 0, 0);
        yset.insert_or_update(0, 0, 0);

        u128 hi = 1;
        u128 tenp = 1;
        u128 ones = 1;
        int leading_start = (start > 0) ? start : n;
        const int rem = n % 10;
        const int scale = rem ? 100 : 10;

        while (true) {
            if (static_cast<u128>(n) % (tenp * 10) != 0) {
                u128 ans = static_cast<u128>(leading_start) * ones;
                for (int f = leading_start; f <= 9; ++f) {
                    const int v = static_cast<int>(ans % static_cast<u128>(n));
                    if (v == 0) {
                        return ans;
                    }

                    bool found = false;
                    u128 best = 0;

                    const int g_begin = rem ? -f : -f;
                    const int g_end = rem ? (10 - f - 1) : (-f);
                    for (int g = g_begin; g <= g_end; ++g) {
                        if (g == 0) {
                            if (found) break;
                            continue;
                        }

                        const GInfo* info = nullptr;
                        for (const GInfo& gi : g_infos) {
                            if (gi.g == g) {
                                info = &gi;
                                break;
                            }
                        }
                        if (info == nullptr) continue;

                        const int rhs = (v == 0) ? 0 : (n - v);
                        if (rhs % info->gcd_ng != 0) continue;

                        const int reduced_rhs = rhs / info->gcd_ng;
                        const int w0 =
                            static_cast<int>((static_cast<long long>(reduced_rhs) * info->inv) %
                                             info->ng);

                        for (int w = w0; w < n; w += info->ng) {
                            for (std::size_t xi = 0; xi < xset.keys.size(); ++xi) {
                                const int x = xset.keys[xi];
                                int y = w - x;
                                if (y < 0) y += n;
                                const int yi = yset.lookup_index(y);
                                if (yi < 0) continue;

                                const u128 sx = (g > 0) ? xset.vals[xi].min_sum : xset.vals[xi].max_sum;
                                const u128 sy = (g > 0)
                                                    ? yset.vals[static_cast<std::size_t>(yi)].min_sum
                                                    : yset.vals[static_cast<std::size_t>(yi)].max_sum;
                                const u128 s = sx + sy;
                                u128 cand = 0;
                                if (g > 0) {
                                    cand = ans + static_cast<u128>(g) * s;
                                } else {
                                    const u128 dec = static_cast<u128>(-g) * s;
                                    if (dec > ans) continue;
                                    cand = ans - dec;
                                }
                                if (!found || cand < best) {
                                    found = true;
                                    best = cand;
                                }
                            }
                        }
                    }

                    if (found) {
                        return best;
                    }
                    ans += ones;
                }

                ResidueSet& target =
                    (yset.size() < static_cast<std::size_t>(scale) * xset.size()) ? yset : xset;
                target.add_shift(static_cast<int>(hi % static_cast<u128>(n)), tenp);
            } else {
                ones = 0;
            }

            leading_start = leading_start / 10;
            if (leading_start == 0) leading_start = 1;
            hi = (hi * 10) % static_cast<u128>(n);
            tenp *= 10;
            ones += tenp;
        }
    }

    std::string d_of_n(int n) const { return to_string_u128(d_u128(n)); }

    u128 D(int k) const {
        if (k <= 0) return 0;

        long hw = sysconf(_SC_NPROCESSORS_ONLN);
        if (hw < 1) hw = 4;
        unsigned thread_count = static_cast<unsigned>(hw);
        if (thread_count > static_cast<unsigned>(k)) thread_count = static_cast<unsigned>(k);
        if (thread_count < 1U) thread_count = 1U;

        if (thread_count == 1U || k < 2000) {
            u128 total = 0;
            for (int n = 1; n <= k; ++n) total += d_u128(n);
            return total;
        }

        std::atomic<int> next(1);
        constexpr int chunk = 64;

        std::vector<pthread_t> tids(static_cast<std::size_t>(thread_count));
        std::vector<WorkerArgs> args(static_cast<std::size_t>(thread_count));
        std::vector<u128> partial(static_cast<std::size_t>(thread_count), 0);

        bool create_failed = false;
        unsigned started = 0;
        for (unsigned t = 0; t < thread_count; ++t) {
            args[static_cast<std::size_t>(t)] =
                WorkerArgs{this, k, &next, chunk, &partial[static_cast<std::size_t>(t)]};
            if (pthread_create(&tids[static_cast<std::size_t>(t)], nullptr, worker_entry,
                               &args[static_cast<std::size_t>(t)]) != 0) {
                create_failed = true;
                break;
            }
            ++started;
        }
        for (unsigned t = 0; t < started; ++t) {
            pthread_join(tids[static_cast<std::size_t>(t)], nullptr);
        }

        if (create_failed) {
            u128 total = 0;
            for (int n = 1; n <= k; ++n) total += d_u128(n);
            return total;
        }

        u128 total = 0;
        for (u128 v : partial) total += v;
        return total;
    }
};

}  // namespace

int main() {
    const Solver solver;

    assert(solver.d_of_n(12) == "12");
    assert(solver.d_of_n(102) == "1122");
    assert(solver.d_of_n(103) == "515");
    assert(solver.d_of_n(290) == "11011010");
    assert(solver.d_of_n(317) == "211122");

    assert(solver.D(110) == 11047);
    assert(solver.D(150) == 53312);
    assert(solver.D(500) == 29570988);

    const u128 answer = solver.D(50'000);
    std::cout << scientific_13sig_u128(answer) << '\n';
    return 0;
}
