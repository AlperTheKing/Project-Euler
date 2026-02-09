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

std::array<u128, 20> build_pow10() {
    std::array<u128, 20> pow10{};
    pow10[0] = 1;
    for (int i = 1; i < static_cast<int>(pow10.size()); ++i) {
        pow10[i] = pow10[i - 1] * static_cast<u128>(10);
    }
    return pow10;
}

u128 g_value(const u64 n, const std::array<u128, 20>& pow10) {
    const std::string s = std::to_string(n);
    const int d = static_cast<int>(s.size());

    std::vector<int> pi(static_cast<std::size_t>(d), 0);
    for (int i = 1; i < d; ++i) {
        int j = pi[static_cast<std::size_t>(i - 1)];
        while (j > 0 && s[static_cast<std::size_t>(i)] != s[static_cast<std::size_t>(j)]) {
            j = pi[static_cast<std::size_t>(j - 1)];
        }
        if (s[static_cast<std::size_t>(i)] == s[static_cast<std::size_t>(j)]) {
            ++j;
        }
        pi[static_cast<std::size_t>(i)] = j;
    }

    u128 expected_end = 0;
    for (int len = d; len > 0; len = pi[static_cast<std::size_t>(len - 1)]) {
        expected_end += pow10[static_cast<std::size_t>(len)];
    }

    return expected_end - static_cast<u128>(d) + 1;
}

u128 solve_sum(const u64 numerator, const u64 max_n, const std::array<u128, 20>& pow10) {
    u128 total = 0;
    for (u64 n = 2; n <= max_n; ++n) {
        total += g_value(numerator / n, pow10);
    }
    return total;
}

bool run_checkpoints() {
    const auto pow10 = build_pow10();

    if (g_value(535ULL, pow10) != static_cast<u128>(1008)) {
        std::cerr << "Checkpoint failed: g(535)\n";
        return false;
    }

    const u128 sample_sum = solve_sum(1000000ULL, 999ULL, pow10);
    if (sample_sum != static_cast<u128>(27280188ULL)) {
        std::cerr << "Checkpoint failed: sample sum\n";
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

    const auto pow10 = build_pow10();
    const u128 answer = solve_sum(10000000000000000ULL, 999999ULL, pow10);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
