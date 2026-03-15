#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 gcd(u64 a, u64 b) {
    while (b != 0) {
        const u64 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

bool is_representable(const u64 n, const u64 a, const u64 b) {
    for (u64 x = 0; x * a <= n; ++x) {
        if ((n - x * a) % b == 0) {
            return true;
        }
    }
    return false;
}

u128 brute_force(const u64 a, const u64 b) {
    std::vector<u64> gaps;
    for (u64 n = 1; n < a * b; ++n) {
        if (!is_representable(n, a, b)) {
            gaps.push_back(n);
        }
    }

    assert(gaps.size() <= 20);

    const int n = static_cast<int>(gaps.size());
    std::vector<std::vector<bool>> attack(static_cast<std::size_t>(n),
                                          std::vector<bool>(static_cast<std::size_t>(n), false));
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const bool bad = is_representable(gaps[static_cast<std::size_t>(j)] -
                                                  gaps[static_cast<std::size_t>(i)],
                                              a,
                                              b);
            attack[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = bad;
            attack[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = bad;
        }
    }

    u128 total = 0;
    const u64 full = 1ULL << n;
    for (u64 mask = 0; mask < full; ++mask) {
        bool ok = true;
        u128 subtotal = 0;
        for (int i = 0; i < n && ok; ++i) {
            if (((mask >> i) & 1ULL) == 0ULL) {
                continue;
            }
            subtotal += static_cast<u128>(gaps[static_cast<std::size_t>(i)]);
            u64 rest = mask & (~0ULL << (i + 1));
            while (rest != 0) {
                const int j = __builtin_ctzll(rest);
                if (attack[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]) {
                    ok = false;
                    break;
                }
                rest &= rest - 1;
            }
        }
        if (ok) {
            total += subtotal;
        }
    }

    return total;
}

u128 solve(u64 a, u64 b) {
    if (a > b) {
        std::swap(a, b);
    }

    const int ai = static_cast<int>(a);
    const int bi = static_cast<int>(b);
    const u64 ab_u64 = a * b;
    const u128 ab = static_cast<u128>(ab_u64);

    std::vector<u128> counts(static_cast<std::size_t>(ai + 1), 0);
    std::vector<u128> sums(static_cast<std::size_t>(ai + 1), 0);
    counts[static_cast<std::size_t>(ai)] = 1;

    for (int x = 1; x < bi; ++x) {
        const int h = static_cast<int>((ab_u64 - 1 - a * static_cast<u64>(x)) / b);
        std::vector<u128> next_counts = counts;
        std::vector<u128> next_sums = sums;

        for (int limit = 1; limit <= ai; ++limit) {
            const u128 count = counts[static_cast<std::size_t>(limit)];
            if (count == 0) {
                continue;
            }

            const u128 prefix_sum = sums[static_cast<std::size_t>(limit)];
            const int ymax = std::min(h, limit - 1);
            for (int y = 1; y <= ymax; ++y) {
                const u128 gap = ab - static_cast<u128>(a) * static_cast<u128>(x) -
                                 static_cast<u128>(b) * static_cast<u128>(y);
                next_counts[static_cast<std::size_t>(y)] += count;
                next_sums[static_cast<std::size_t>(y)] += prefix_sum + count * gap;
            }
        }

        counts.swap(next_counts);
        sums.swap(next_sums);
    }

    u128 answer = 0;
    for (const u128 value : sums) {
        answer += value;
    }
    return answer;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

void run_checkpoints() {
    for (u64 a = 2; a <= 7; ++a) {
        for (u64 b = a + 1; b <= 8; ++b) {
            if (gcd(a, b) != 1 || a * b > 35) {
                continue;
            }
            assert(solve(a, b) == brute_force(a, b));
            assert(solve(a, b) == solve(b, a));
        }
    }

    assert(solve(3, 5) == 23);
    assert(solve(5, 13) == 16336);
}

}  // namespace

int main(int argc, char** argv) {
    bool should_run_checkpoints = true;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--skip-checkpoints") {
            should_run_checkpoints = false;
        }
    }

    if (should_run_checkpoints) {
        run_checkpoints();
    }

    std::cout << to_string_u128(solve(19, 53)) << '\n';
    return 0;
}
