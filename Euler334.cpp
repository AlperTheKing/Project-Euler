#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i128 = __int128_t;
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

std::vector<u64> generate_beans(const int count) {
    std::vector<u64> beans;
    beans.reserve(static_cast<std::size_t>(count));

    u64 t = 123456ULL;
    for (int i = 0; i < count; ++i) {
        if ((t & 1ULL) == 0ULL) {
            t /= 2ULL;
        } else {
            t = (t / 2ULL) ^ 926252ULL;
        }
        const u64 b = (t & ((1ULL << 11) - 1ULL)) + 1ULL;
        beans.push_back(b);
    }

    return beans;
}

i128 floor_div(const i128 a, const i128 b) {
    i128 q = a / b;
    i128 r = a % b;
    if (r < 0) {
        q -= 1;
    }
    return q;
}

i128 sum_squares_arith(const i128 first, const i128 n) {
    // Sum_{k=0}^{n-1} (first + k)^2
    if (n <= 0) {
        return 0;
    }
    const i128 sum_k = n * (n - 1) / 2;
    const i128 sum_k2 = n * (n - 1) * (2 * n - 1) / 6;
    return n * first * first + 2 * first * sum_k + sum_k2;
}

u128 moves_needed(const std::vector<u64>& beans) {
    i128 total = 0;
    i128 first_moment = 0;
    i128 second_moment = 0;

    for (i128 i = 0; i < static_cast<i128>(beans.size()); ++i) {
        const i128 b = static_cast<i128>(beans[static_cast<std::size_t>(i)]);
        total += b;
        first_moment += i * b;
        second_moment += i * i * b;
    }

    // In every move, second moment increases by exactly 2.
    // Stable final states are 0/1-occupancy sets with same total and first moment;
    // the attained stabilization is the minimum-second-moment such set.
    const i128 B = total;
    const i128 target_shift = first_moment - B * (B - 1) / 2;

    const i128 q = floor_div(target_shift, B);
    const i128 r = target_shift - q * B;  // 0 <= r < B

    const i128 n1 = B - r;
    const i128 n2 = r;

    // Final occupied positions:
    // q, q+1, ..., q+n1-1, q+n1+1, ..., q+B
    const i128 q_final = sum_squares_arith(q, n1) + sum_squares_arith(q + n1 + 1, n2);

    const i128 diff = q_final - second_moment;
    return static_cast<u128>(diff / 2);
}

bool run_checkpoints() {
    if (moves_needed(std::vector<u64>{2ULL, 3ULL}) != static_cast<u128>(8ULL)) {
        std::cerr << "Checkpoint failed: [2,3] should take 8 moves\n";
        return false;
    }

    const std::vector<u64> b2 = generate_beans(2);
    if (b2.size() != 2 || b2[0] != 289ULL || b2[1] != 145ULL) {
        std::cerr << "Checkpoint failed: b1/b2 generation\n";
        return false;
    }

    if (moves_needed(b2) != static_cast<u128>(3419100ULL)) {
        std::cerr << "Checkpoint failed: two-bowl sample moves\n";
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

    const std::vector<u64> beans = generate_beans(1500);
    const u128 answer = moves_needed(beans);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
