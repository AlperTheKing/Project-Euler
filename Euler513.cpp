#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

struct Options {
    int n = 100000;
    i64 segment_size = 1000000;
    unsigned requested_threads = 0U;
    bool run_checkpoints = true;
};

struct CEntry {
    u32 key;
    u32 c;
};

struct ABEntry {
    u32 key;
    u32 b;
    u32 sum;
};

bool parse_int_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(ch - '0');
        if (parsed > static_cast<i64>(std::numeric_limits<int>::max())) {
            return false;
        }
    }

    value = static_cast<int>(parsed);
    return true;
}

bool parse_i64_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const int digit = ch - '0';
        if (parsed > (std::numeric_limits<i64>::max() - digit) / 10) {
            return false;
        }
        parsed = parsed * 10 + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const std::string& prefix,
                                 unsigned& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(ch - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
        if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
            return false;
        }
    }

    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        if (parse_i64_after_prefix(arg, "--segment=", options.segment_size)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.n < 1) {
        std::cerr << "--n must be >= 1\n";
        return false;
    }
    if (options.segment_size < 1) {
        std::cerr << "--segment must be >= 1\n";
        return false;
    }

    return true;
}

unsigned pick_thread_count(const unsigned requested) {
    if (requested > 0U) {
        return requested;
    }
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0U) {
        hw = 4U;
    }
    return hw;
}

i64 isqrt_floor(const i64 x) {
    i64 r = static_cast<i64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r > x) {
        --r;
    }
    return r;
}

struct Solver {
    int n = 1;
    int m = 0;
    i64 max_n = 0;
    i64 segment_size = 1;

    u64 solve_segment(const i64 L,
                      const i64 R,
                      std::vector<CEntry>& c_entries,
                      std::vector<ABEntry>& ab_entries) const {
        c_entries.clear();
        ab_entries.clear();

        for (int x = 1; x <= m; ++x) {
            const i64 xx = 1LL * x * x;
            if (xx + xx > R) {
                break;
            }

            const i64 lo = (L + 1) - xx;
            int t_min = x;
            if (lo > 0) {
                i64 r = isqrt_floor(lo);
                if (r * r < lo) {
                    ++r;
                }
                if (r > t_min) {
                    t_min = static_cast<int>(r);
                }
            }

            const i64 hi = R - xx;
            if (hi < 1LL * t_min * t_min) {
                continue;
            }
            const int t_max = static_cast<int>(isqrt_floor(hi));
            if (t_max < t_min) {
                continue;
            }

            for (int t = t_min; t <= t_max; ++t) {
                const i64 n_value = xx + 1LL * t * t;
                const u32 key = static_cast<u32>(n_value - L);
                c_entries.push_back({key, static_cast<u32>(x)});
                if (t != x && t <= m) {
                    c_entries.push_back({key, static_cast<u32>(t)});
                }
            }
        }

        const i64 sum_lo = 2LL * L + 2;
        const i64 sum_hi = 2LL * R;
        for (int a = 1; a <= n; ++a) {
            const i64 aa = 1LL * a * a;
            if (aa + aa > sum_hi) {
                break;
            }

            const i64 lo = sum_lo - aa;
            int b_min = a;
            if (lo > 0) {
                i64 r = isqrt_floor(lo);
                if (r * r < lo) {
                    ++r;
                }
                if (r > b_min) {
                    b_min = static_cast<int>(r);
                }
            }
            if ((b_min & 1) != (a & 1)) {
                ++b_min;
            }

            const i64 hi = sum_hi - aa;
            if (hi < 1LL * b_min * b_min) {
                continue;
            }
            int b_max = static_cast<int>(isqrt_floor(hi));
            if (b_max > n) {
                b_max = n;
            }
            if (b_max < b_min) {
                continue;
            }

            for (int b = b_min; b <= b_max; b += 2) {
                const i64 s = aa + 1LL * b * b;
                if (s & 1LL) {
                    continue;
                }
                const i64 n_value = s / 2LL;
                if (n_value <= L || n_value > R) {
                    continue;
                }
                const u32 key = static_cast<u32>(n_value - L);
                ab_entries.push_back(
                    {key, static_cast<u32>(b), static_cast<u32>(a + b)});
            }
        }

        std::sort(c_entries.begin(), c_entries.end(),
                  [](const CEntry& lhs, const CEntry& rhs) {
                      if (lhs.key != rhs.key) {
                          return lhs.key < rhs.key;
                      }
                      return lhs.c < rhs.c;
                  });

        std::sort(ab_entries.begin(), ab_entries.end(),
                  [](const ABEntry& lhs, const ABEntry& rhs) {
                      if (lhs.key != rhs.key) {
                          return lhs.key < rhs.key;
                      }
                      if (lhs.b != rhs.b) {
                          return lhs.b < rhs.b;
                      }
                      return lhs.sum < rhs.sum;
                  });

        u64 segment_answer = 0ULL;
        std::size_t i = 0U;
        std::size_t j = 0U;

        while (i < c_entries.size() && j < ab_entries.size()) {
            if (c_entries[i].key < ab_entries[j].key) {
                const u32 key = c_entries[i].key;
                while (i < c_entries.size() && c_entries[i].key == key) {
                    ++i;
                }
                continue;
            }
            if (ab_entries[j].key < c_entries[i].key) {
                const u32 key = ab_entries[j].key;
                while (j < ab_entries.size() && ab_entries[j].key == key) {
                    ++j;
                }
                continue;
            }

            const u32 key = c_entries[i].key;
            const std::size_t i_begin = i;
            while (i < c_entries.size() && c_entries[i].key == key) {
                ++i;
            }
            const std::size_t i_end = i;

            const std::size_t j_begin = j;
            while (j < ab_entries.size() && ab_entries[j].key == key) {
                ++j;
            }
            const std::size_t j_end = j;

            for (std::size_t idx = j_begin; idx < j_end; ++idx) {
                const int lo_c = (static_cast<int>(ab_entries[idx].b) + 1) / 2;
                const int hi_c = (static_cast<int>(ab_entries[idx].sum) - 1) / 2;
                if (lo_c > hi_c) {
                    continue;
                }

                const auto lower = std::lower_bound(
                    c_entries.begin() + static_cast<std::ptrdiff_t>(i_begin),
                    c_entries.begin() + static_cast<std::ptrdiff_t>(i_end),
                    static_cast<u32>(lo_c),
                    [](const CEntry& entry, const u32 value) {
                        return entry.c < value;
                    });
                const auto upper = std::upper_bound(
                    c_entries.begin() + static_cast<std::ptrdiff_t>(i_begin),
                    c_entries.begin() + static_cast<std::ptrdiff_t>(i_end),
                    static_cast<u32>(hi_c),
                    [](const u32 value, const CEntry& entry) {
                        return value < entry.c;
                    });

                segment_answer += static_cast<u64>(upper - lower);
            }
        }

        return segment_answer;
    }

    u64 solve(const unsigned thread_count) const {
        const i64 segment_count = (max_n + segment_size - 1) / segment_size;
        std::atomic<i64> next_segment(0);

        std::vector<u64> partial(thread_count, 0ULL);
        std::vector<std::thread> pool;
        pool.reserve(thread_count);

        for (unsigned t = 0U; t < thread_count; ++t) {
            pool.emplace_back([&, t]() {
                std::vector<CEntry> c_entries;
                std::vector<ABEntry> ab_entries;
                c_entries.reserve(1U << 20U);
                ab_entries.reserve(1U << 20U);

                u64 local = 0ULL;
                while (true) {
                    const i64 idx = next_segment.fetch_add(1, std::memory_order_relaxed);
                    if (idx >= segment_count) {
                        break;
                    }

                    const i64 L = idx * segment_size;
                    const i64 R = std::min(max_n, L + segment_size);
                    local += solve_segment(L, R, c_entries, ab_entries);
                }

                partial[t] = local;
            });
        }

        for (std::thread& th : pool) {
            th.join();
        }

        u64 answer = 0ULL;
        for (u64 value : partial) {
            answer += value;
        }
        return answer;
    }
};

u64 compute_f(const int n, const i64 segment_size, const unsigned thread_count) {
    const int m = n / 2;
    Solver solver{n, m, 4LL * m * m, segment_size};
    return solver.solve(thread_count);
}

bool run_checkpoints(const i64 segment_size, const unsigned thread_count) {
    struct Checkpoint {
        int n;
        u64 expected;
    };

    const Checkpoint checkpoints[] = {
        {10, 3ULL},
        {50, 165ULL},
    };

    for (const Checkpoint& checkpoint : checkpoints) {
        const u64 got = compute_f(checkpoint.n, segment_size, thread_count);
        if (got != checkpoint.expected) {
            std::cerr << "Validation failed for n=" << checkpoint.n
                      << ": got " << got
                      << ", expected " << checkpoint.expected << '\n';
            return false;
        }
    }

    std::cerr << "Validation checkpoints passed.\n";
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const unsigned threads = pick_thread_count(options.requested_threads);
    if (options.run_checkpoints &&
        !run_checkpoints(options.segment_size, std::min(threads, 4U))) {
        return 1;
    }

    std::cout << compute_f(options.n, options.segment_size, threads) << '\n';
    return 0;
}
