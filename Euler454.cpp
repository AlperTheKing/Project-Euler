#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = __uint128_t;

struct Options {
    u64 l = 1'000'000'000'000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = static_cast<u64>(std::stoull(tail));
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--l=", options.l)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.l >= 1ULL;
}

u64 integer_sqrt(const u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        const unsigned digit = static_cast<unsigned>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

u64 solve(const u64 L) {
    const u64 root = integer_sqrt(1ULL + 4ULL * L);
    const u32 max_b = static_cast<u32>((root - 1ULL) / 2ULL);

    std::vector<int> mu(static_cast<std::size_t>(max_b) + 1U, 0);
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(max_b / 10U));
    std::vector<bool> is_composite(static_cast<std::size_t>(max_b) + 1U, false);
    mu[1] = 1;

    for (u32 i = 2U; i <= max_b; ++i) {
        if (!is_composite[i]) {
            primes.push_back(i);
            mu[i] = -1;
        }
        for (const u32 p : primes) {
            const u64 v = static_cast<u64>(i) * static_cast<u64>(p);
            if (v > max_b) {
                break;
            }
            is_composite[static_cast<std::size_t>(v)] = true;
            if (i % p == 0U) {
                mu[static_cast<std::size_t>(v)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(v)] = -mu[i];
        }
    }

    std::vector<u32> divisor_count(static_cast<std::size_t>(max_b) + 1U, 0U);
    for (u32 d = 1U; d <= max_b; ++d) {
        if (mu[d] == 0) {
            continue;
        }
        for (u32 b = d; b <= max_b; b += d) {
            ++divisor_count[b];
        }
    }

    std::vector<u32> offset(static_cast<std::size_t>(max_b) + 2U, 0U);
    for (u32 b = 1U; b <= max_b; ++b) {
        offset[static_cast<std::size_t>(b) + 1U] =
            offset[static_cast<std::size_t>(b)] + divisor_count[b];
    }

    std::vector<int> signed_divisor(offset[static_cast<std::size_t>(max_b) + 1U], 0);
    std::vector<u32> fill = offset;

    for (u32 d = 1U; d <= max_b; ++d) {
        const int mud = mu[d];
        if (mud == 0) {
            continue;
        }
        const int sd = mud * static_cast<int>(d);
        for (u32 b = d; b <= max_b; b += d) {
            signed_divisor[fill[b]++] = sd;
        }
    }

    u128 answer = 0;
    for (u32 b = 1U; b <= max_b; ++b) {
        const u64 M = L / b;
        const u64 upper = std::min<u64>(2ULL * b - 1ULL, M);
        if (upper <= b) {
            continue;
        }

        u64 c = static_cast<u64>(b) + 1ULL;
        while (c <= upper) {
            const u64 q = M / c;
            const u64 c2 = std::min<u64>(upper, M / q);

            const u64 left_minus_1 = c - 1ULL;
            i64 coprime_count = 0;
            for (u32 idx = offset[b]; idx < offset[static_cast<std::size_t>(b) + 1U]; ++idx) {
                const int sd = signed_divisor[idx];
                const u32 d = static_cast<u32>(sd > 0 ? sd : -sd);
                const i64 term =
                    static_cast<i64>(c2 / d) - static_cast<i64>(left_minus_1 / d);
                coprime_count += (sd > 0) ? term : -term;
            }

            answer += static_cast<u128>(q) * static_cast<u128>(coprime_count);
            c = c2 + 1ULL;
        }
    }

    return static_cast<u64>(answer);
}

bool run_checkpoints() {
    if (solve(15ULL) != 4ULL) {
        std::cerr << "Checkpoint failed: L(15)\n";
        return false;
    }
    if (solve(1000ULL) != 1069ULL) {
        std::cerr << "Checkpoint failed: L(1000)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const u64 answer = solve(options.l);
    std::cout << answer << '\n';
    return 0;
}
