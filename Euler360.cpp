#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u8 = std::uint8_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 isqrt_u64(const u64 n) {
    if (n == 0ULL) {
        return 0ULL;
    }
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
}

u64 pow_u64(u64 base, unsigned exp) {
    u64 result = 1;
    while (exp > 0U) {
        if ((exp & 1U) != 0U) {
            result *= base;
        }
        base *= base;
        exp >>= 1U;
    }
    return result;
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }
    std::string out;
    while (value > 0U) {
        const unsigned digit = static_cast<unsigned>(value % 10U);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

std::vector<u32> build_spf(const u32 limit) {
    std::vector<u32> spf(static_cast<std::size_t>(limit + 1U));
    for (u32 i = 0; i <= limit; ++i) {
        spf[static_cast<std::size_t>(i)] = i;
    }

    for (u32 i = 2; static_cast<u64>(i) * i <= limit; ++i) {
        if (spf[static_cast<std::size_t>(i)] != i) {
            continue;
        }
        for (u32 j = i * i; j <= limit; j += i) {
            if (spf[static_cast<std::size_t>(j)] == j) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }

    return spf;
}

struct FactorTables {
    std::vector<u8> bad_3mod4_odd;
    std::vector<u8> exp5;
    std::vector<u32> one_mod4_factor;
};

FactorTables build_factor_tables(const u32 limit, const std::vector<u32>& spf) {
    FactorTables data;
    data.bad_3mod4_odd.assign(static_cast<std::size_t>(limit + 1U), 0U);
    data.exp5.assign(static_cast<std::size_t>(limit + 1U), 0U);
    data.one_mod4_factor.assign(static_cast<std::size_t>(limit + 1U), 1U);

    for (u32 n = 1; n <= limit; ++n) {
        u32 x = n;
        bool bad = false;
        u8 e5 = 0;
        u32 mult = 1;

        while (x > 1U) {
            const u32 p = spf[static_cast<std::size_t>(x)];
            u32 e = 0;
            do {
                x /= p;
                ++e;
            } while (x % p == 0U);

            if (p == 5U) {
                e5 = static_cast<u8>(e);
            } else if ((p & 3U) == 3U) {
                if ((e & 1U) != 0U) {
                    bad = true;
                }
            } else if ((p & 3U) == 1U) {
                mult *= (e + 1U);
            }
        }

        data.bad_3mod4_odd[static_cast<std::size_t>(n)] = bad ? 1U : 0U;
        data.exp5[static_cast<std::size_t>(n)] = e5;
        data.one_mod4_factor[static_cast<std::size_t>(n)] = mult;
    }

    return data;
}

u128 sum_for_radius_power_of_five(const unsigned exponent) {
    const u64 r = pow_u64(5ULL, exponent);
    const u32 limit = static_cast<u32>(2ULL * r);

    const std::vector<u32> spf = build_spf(limit);
    const FactorTables tables = build_factor_tables(limit, spf);

    const auto& bad = tables.bad_3mod4_odd;
    const auto& e5 = tables.exp5;
    const auto& f1 = tables.one_mod4_factor;

    const u64 work_items = (r >= 1ULL) ? (r - 1ULL) : 0ULL;
    unsigned thread_count = 1U;
    if (work_items >= 1000000ULL) {
        thread_count = std::thread::hardware_concurrency();
        if (thread_count == 0U) {
            thread_count = 4U;
        }
        thread_count = static_cast<unsigned>(std::min<u64>(thread_count, work_items));
    }

    std::vector<u128> partial(thread_count, 0U);
    std::vector<std::thread> workers;
    workers.reserve(thread_count);

    const u64 chunk = (work_items + thread_count - 1ULL) / thread_count;
    for (unsigned t = 0; t < thread_count; ++t) {
        const u64 begin = 1ULL + static_cast<u64>(t) * chunk;
        const u64 end = std::min<u64>(r, begin + chunk);

        workers.emplace_back([&, t, begin, end]() {
            u128 local = 0;
            for (u64 x = begin; x < end; ++x) {
                const u32 a = static_cast<u32>(r - x);
                const u32 b = static_cast<u32>(r + x);

                if (bad[static_cast<std::size_t>(a)] != 0U || bad[static_cast<std::size_t>(b)] != 0U) {
                    continue;
                }

                const u128 ways = 4U * static_cast<u128>(f1[static_cast<std::size_t>(a)]) *
                                  static_cast<u128>(f1[static_cast<std::size_t>(b)]) *
                                  static_cast<u128>(static_cast<u32>(e5[static_cast<std::size_t>(a)]) +
                                                    static_cast<u32>(e5[static_cast<std::size_t>(b)]) + 1U);

                local += 6U * static_cast<u128>(x) * ways;
            }
            partial[t] = local;
        });
    }

    for (auto& th : workers) {
        th.join();
    }

    // x = r contributes r_2(0) = 1.
    u128 total = 6U * static_cast<u128>(r);
    for (const u128 v : partial) {
        total += v;
    }
    return total;
}

u128 sum_for_radius_2a_5b(const unsigned exp2, const unsigned exp5) {
    u128 total = sum_for_radius_power_of_five(exp5);
    total <<= exp2;
    return total;
}

u64 brute_sum(const int r) {
    const int rr = r * r;
    u128 total = 0;

    for (int x = -r; x <= r; ++x) {
        for (int y = -r; y <= r; ++y) {
            const int z2 = rr - x * x - y * y;
            if (z2 < 0) {
                continue;
            }

            const int z = static_cast<int>(isqrt_u64(static_cast<u64>(z2)));
            if (z * z != z2) {
                continue;
            }

            const u64 base = static_cast<u64>(std::abs(x)) + static_cast<u64>(std::abs(y));
            if (z == 0) {
                total += base;
            } else {
                total += base + static_cast<u64>(z);
                total += base + static_cast<u64>(z);
            }
        }
    }

    return static_cast<u64>(total);
}

bool run_checkpoints() {
    if (brute_sum(45) != 34518ULL) {
        std::cerr << "Checkpoint failed: S(45)\n";
        return false;
    }

    if (sum_for_radius_2a_5b(0U, 1U) != 198U) {
        std::cerr << "Checkpoint failed: S(5)\n";
        return false;
    }

    if (sum_for_radius_2a_5b(1U, 1U) != 396U) {
        std::cerr << "Checkpoint failed: S(10)\n";
        return false;
    }

    if (sum_for_radius_2a_5b(0U, 2U) != 5766U) {
        std::cerr << "Checkpoint failed: S(25)\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    // 10^10 = 2^10 * 5^10
    const u128 answer = sum_for_radius_2a_5b(10U, 10U);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
