#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

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

std::vector<u32> precompute_golomb_until(const u64 max_query_position) {
    // 1-indexed: golomb[n] = G(n)
    std::vector<u32> golomb;
    golomb.reserve(11000000);
    golomb.push_back(0U);
    golomb.push_back(1U);

    u128 products = 1;  // sum_{i=1..k} i * G(i)

    for (u64 i = 2; products < static_cast<u128>(max_query_position); ++i) {
        const u64 prev = i - 1;
        const u32 g_prev = golomb[static_cast<std::size_t>(prev)];
        const u32 nested = golomb[static_cast<std::size_t>(g_prev)];
        const u64 idx = i - static_cast<u64>(nested);
        const u32 current = static_cast<u32>(1U + golomb[static_cast<std::size_t>(idx)]);

        golomb.push_back(current);
        products += static_cast<u128>(i) * static_cast<u128>(current);
    }

    return golomb;
}

struct QueryState {
    u64 index = 1;      // current k
    u128 sum = 1;       // S(k) = sum G(i)
    u128 products = 1;  // P(k) = sum i*G(i)

    u128 prev_sum = 0;       // S(k-1)
    u128 prev_products = 0;  // P(k-1)
};

u64 golomb_at(const u64 x, const std::vector<u32>& golomb, QueryState& st) {
    while (st.products < static_cast<u128>(x)) {
        ++st.index;
        st.prev_sum = st.sum;
        st.prev_products = st.products;

        const u32 g = golomb[static_cast<std::size_t>(st.index)];
        st.sum += static_cast<u128>(g);
        st.products += static_cast<u128>(st.index) * static_cast<u128>(g);
    }

    const u128 span_index = static_cast<u128>(st.index);
    const u128 offset = (static_cast<u128>(x) - st.prev_products + span_index - 1U) / span_index;
    return static_cast<u64>(st.prev_sum + offset);
}

u128 sum_golomb_cubes(const u64 limit, const std::vector<u32>& golomb) {
    QueryState st;
    u128 sum = 0;

    for (u64 n = 1; n < limit; ++n) {
        const u64 cube = n * n * n;
        sum += static_cast<u128>(golomb_at(cube, golomb, st));
    }

    return sum;
}

bool run_checkpoints(const std::vector<u32>& golomb) {
    {
        QueryState st;
        if (golomb_at(1000ULL, golomb, st) != 86ULL) {
            std::cerr << "Checkpoint failed: G(10^3)\n";
            return false;
        }
    }

    {
        QueryState st;
        if (golomb_at(1000000ULL, golomb, st) != 6137ULL) {
            std::cerr << "Checkpoint failed: G(10^6)\n";
            return false;
        }
    }

    {
        const u128 sample = sum_golomb_cubes(1000ULL, golomb);
        if (sample != static_cast<u128>(153506976ULL)) {
            std::cerr << "Checkpoint failed: sum_{1<=n<10^3} G(n^3)\n";
            return false;
        }
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

    const u64 limit = 1000000ULL;
    const u64 max_query_position = (limit - 1ULL) * (limit - 1ULL) * (limit - 1ULL);
    const std::vector<u32> golomb = precompute_golomb_until(max_query_position);

    if (!skip_checkpoints && !run_checkpoints(golomb)) {
        return 2;
    }

    const u128 answer = sum_golomb_cubes(limit, golomb);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
