#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    bool run_checkpoints = true;
};

const std::string A100 =
    "14159265358979323846264338327950288419716939937510"
    "58209749445923078164062862089986280348253421170679";
const std::string B100 =
    "82148086513282306647093844609550582231725359408128"
    "48111745028410270193852110555964462294895493038196";

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

int digit_at(const std::string& A, const std::string& B, u64 n) {
    std::vector<u64> len = {0ULL, static_cast<u64>(A.size()), static_cast<u64>(B.size())};
    while (len.back() < n) {
        const int m = static_cast<int>(len.size()) - 1;
        len.push_back(len[static_cast<std::size_t>(m - 1)] + len[static_cast<std::size_t>(m)]);
    }

    int k = static_cast<int>(len.size()) - 1;
    while (k > 2) {
        if (n <= len[static_cast<std::size_t>(k - 2)]) {
            k -= 2;
        } else {
            n -= len[static_cast<std::size_t>(k - 2)];
            k -= 1;
        }
    }

    return (k == 1) ? (A[static_cast<std::size_t>(n - 1)] - '0')
                    : (B[static_cast<std::size_t>(n - 1)] - '0');
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

u128 solve() {
    u128 answer = 0;
    u128 pow10 = 1;

    for (int n = 0; n <= 17; ++n) {
        u64 index = static_cast<u64>(127 + 19 * n);
        for (int i = 0; i < n; ++i) {
            index *= 7ULL;
        }
        const int digit = digit_at(A100, B100, index);
        answer += pow10 * static_cast<u128>(digit);
        pow10 *= 10U;
    }

    return answer;
}

bool run_checkpoints() {
    const std::string A = "1415926535";
    const std::string B = "8979323846";
    if (digit_at(A, B, 35ULL) != 9) {
        std::cerr << "Checkpoint failed for statement example" << '\n';
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

    std::cout << to_string_u128(solve()) << '\n';
    return 0;
}
