#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 50515093ULL;
constexpr u64 kTargetN = 2000000000ULL;

std::vector<u32> build_period_sequence() {
    u64 s = 290797ULL;
    s = (s * s) % kMod;  // S_1
    const u32 first = static_cast<u32>(s);

    std::vector<u32> seq;
    seq.reserve(6500000);
    while (true) {
        seq.push_back(static_cast<u32>(s));
        s = (s * s) % kMod;
        if (static_cast<u32>(s) == first) {
            break;
        }
    }
    return seq;
}

void build_spans(const std::vector<u32>& a, std::vector<u32>& d_left, std::vector<u32>& d_right) {
    const u32 p = static_cast<u32>(a.size());
    d_left.assign(static_cast<std::size_t>(p), 0U);   // 0 means no strictly smaller in previous period
    d_right.assign(static_cast<std::size_t>(p), 0U);  // always positive

    auto value_at = [&](const u32 pos) -> u32 {
        return a[static_cast<std::size_t>((pos - 1U) % p)];
    };

    std::vector<u32> stack;
    stack.reserve(static_cast<std::size_t>(2U * p + 8U));

    // Next position with value <= current.
    for (u32 pos = 2U * p; pos >= 1U; --pos) {
        const u32 v = value_at(pos);
        while (!stack.empty() && value_at(stack.back()) > v) {
            stack.pop_back();
        }
        if (pos <= p) {
            d_right[static_cast<std::size_t>(pos - 1U)] = stack.back() - pos;
        }
        stack.push_back(pos);
        if (pos == 1U) {
            break;
        }
    }

    // Previous position with value < current, limited to at most one period back.
    stack.clear();
    for (u32 pos = 1U; pos <= 2U * p; ++pos) {
        const u32 v = value_at(pos);
        while (!stack.empty() && value_at(stack.back()) >= v) {
            stack.pop_back();
        }
        if (pos > p) {
            const u32 idx = pos - p;  // 1-based index in one period
            if (!stack.empty() && stack.back() >= idx) {
                d_left[static_cast<std::size_t>(idx - 1U)] = pos - stack.back();
            } else {
                d_left[static_cast<std::size_t>(idx - 1U)] = 0U;
            }
        }
        stack.push_back(pos);
    }
}

u128 sum_subarray_min_prefix(const std::vector<u32>& a,
                             const std::vector<u32>& d_left,
                             const std::vector<u32>& d_right,
                             const u64 n) {
    const u64 p = static_cast<u64>(a.size());
    const u64 q = n / p;
    const u64 r = n % p;

    u128 total = 0U;
    for (u64 idx = 1ULL; idx <= p; ++idx) {
        const u64 occ = q + ((idx <= r) ? 1ULL : 0ULL);
        if (occ == 0ULL) {
            continue;
        }

        const u64 v = static_cast<u64>(a[static_cast<std::size_t>(idx - 1ULL)]);
        const u64 dr = static_cast<u64>(d_right[static_cast<std::size_t>(idx - 1ULL)]);
        const u64 dl = static_cast<u64>(d_left[static_cast<std::size_t>(idx - 1ULL)]);

        const u64 rem_last = (idx <= r) ? (r - idx + 1ULL) : (p + r - idx + 1ULL);
        const u64 r_last = std::min(dr, rem_last);

        if (dl == 0ULL) {
            // No strictly smaller element exists to the left in any previous period.
            if (occ == 1ULL) {
                total += static_cast<u128>(v) * static_cast<u128>(idx) * static_cast<u128>(r_last);
            } else {
                const u64 n1 = occ - 1ULL;
                const u128 sum_l = static_cast<u128>(n1) * static_cast<u128>(idx) +
                                   static_cast<u128>(p) * static_cast<u128>(n1 - 1ULL) * static_cast<u128>(n1) / 2U;
                total += static_cast<u128>(v) * static_cast<u128>(dr) * sum_l;

                const u64 l_last = (occ - 1ULL) * p + idx;
                total += static_cast<u128>(v) * static_cast<u128>(l_last) * static_cast<u128>(r_last);
            }
        } else {
            if (occ == 1ULL) {
                const u64 l0 = std::min(dl, idx);
                total += static_cast<u128>(v) * static_cast<u128>(l0) * static_cast<u128>(r_last);
            } else {
                const u64 l0 = std::min(dl, idx);
                total += static_cast<u128>(v) * static_cast<u128>(l0) * static_cast<u128>(dr);
                if (occ > 2ULL) {
                    total += static_cast<u128>(v) * static_cast<u128>(dl) * static_cast<u128>(dr) *
                             static_cast<u128>(occ - 2ULL);
                }
                total += static_cast<u128>(v) * static_cast<u128>(dl) * static_cast<u128>(r_last);
            }
        }
    }

    return total;
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

bool run_checkpoints(const std::vector<u32>& a, const std::vector<u32>& d_left, const std::vector<u32>& d_right) {
    if (a.size() != 6308948U) {
        std::cerr << "Checkpoint failed: period length\n";
        return false;
    }

    if (sum_subarray_min_prefix(a, d_left, d_right, 10ULL) != 432256955ULL) {
        std::cerr << "Checkpoint failed: M(10)\n";
        return false;
    }

    if (sum_subarray_min_prefix(a, d_left, d_right, 10000ULL) != 3264567774119ULL) {
        std::cerr << "Checkpoint failed: M(10000)\n";
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

    const std::vector<u32> period = build_period_sequence();
    std::vector<u32> d_left;
    std::vector<u32> d_right;
    build_spans(period, d_left, d_right);

    if (!skip_checkpoints && !run_checkpoints(period, d_left, d_right)) {
        return 2;
    }

    const u128 answer = sum_subarray_min_prefix(period, d_left, d_right, kTargetN);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
