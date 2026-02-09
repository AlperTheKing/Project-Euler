#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

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

u128 count_pairs(const u64 N, const int M) {
    const u64 period = 6ULL * static_cast<u64>(M);
    const u64 full_cycles = N / period;
    const u64 remainder = N % period;

    std::vector<u64> freq_in_period(static_cast<std::size_t>(M), 0ULL);
    std::vector<u64> freq_in_tail(static_cast<std::size_t>(M), 0ULL);

    u128 weighted_sum = 0;  // sum_{k=1}^{n} k * a_k
    int prefix_mod = 0;     // sum_{i=1}^{n} a_i (mod M)

    for (u64 n = 0; n < period; ++n) {
        ++freq_in_period[static_cast<std::size_t>(prefix_mod)];
        if (n <= remainder) {
            ++freq_in_tail[static_cast<std::size_t>(prefix_mod)];
        }

        const u64 idx = n + 1ULL;
        const u64 a = (idx == 1ULL) ? 1ULL : static_cast<u64>(weighted_sum % idx);
        weighted_sum += static_cast<u128>(idx) * static_cast<u128>(a);

        prefix_mod += static_cast<int>(a % static_cast<u64>(M));
        if (prefix_mod >= M) {
            prefix_mod %= M;
        }
    }

    u128 answer = 0;
    for (int r = 0; r < M; ++r) {
        const u128 count = static_cast<u128>(full_cycles) * static_cast<u128>(freq_in_period[static_cast<std::size_t>(r)]) +
                           static_cast<u128>(freq_in_tail[static_cast<std::size_t>(r)]);
        answer += count * (count - 1) / 2;
    }

    return answer;
}

bool run_checkpoints() {
    {
        // Statement: first 10 values of a_n.
        std::array<u64, 10> expected = {1, 1, 0, 3, 0, 3, 5, 4, 1, 9};
        std::array<u64, 10> got{};

        u128 weighted_sum = 0;
        for (u64 n = 1; n <= 10; ++n) {
            const u64 a = (n == 1ULL) ? 1ULL : static_cast<u64>(weighted_sum % n);
            got[static_cast<std::size_t>(n - 1)] = a;
            weighted_sum += static_cast<u128>(n) * static_cast<u128>(a);
        }

        if (got != expected) {
            std::cerr << "Checkpoint failed: first 10 sequence elements mismatch\n";
            return false;
        }
    }

    if (count_pairs(10ULL, 10) != static_cast<u128>(4ULL)) {
        std::cerr << "Checkpoint failed: f(10,10)\n";
        return false;
    }

    if (count_pairs(10000ULL, 1000) != static_cast<u128>(97158ULL)) {
        std::cerr << "Checkpoint failed: f(10^4,10^3)\n";
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

    const u128 answer = count_pairs(1000000000000ULL, 1000000);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
