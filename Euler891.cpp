#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

namespace {

constexpr std::int64_t kCycleSeconds = 43200;

std::int64_t extended_gcd(std::int64_t a, std::int64_t b, std::int64_t& x,
                          std::int64_t& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    std::int64_t x1 = 0;
    std::int64_t y1 = 0;
    std::int64_t g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

}  // namespace

int main() {
    const std::int64_t speeds[3] = {1, 12, 720};
    std::array<int, 3> idx = {0, 1, 2};

    std::vector<std::pair<std::int64_t, std::int64_t>> moments;
    moments.reserve(2100000);

    // For each non-identity permutation, eliminate rotation via differences to
    // get a linear congruence system: t = 43200*n/|det| and n' = C*n (mod |det|).
    do {
        if (idx[0] == 0 && idx[1] == 1 && idx[2] == 2) {
            continue;
        }
        std::int64_t pa = speeds[idx[0]];
        std::int64_t pb = speeds[idx[1]];
        std::int64_t pc = speeds[idx[2]];

        std::int64_t A = pb - pa;
        std::int64_t B = pc - pa;
        std::int64_t det = 719 * A - 11 * B;
        if (det == 0) {
            continue;
        }
        std::int64_t D = det > 0 ? det : -det;

        std::int64_t absA = A >= 0 ? A : -A;
        std::int64_t absB = B >= 0 ? B : -B;
        std::int64_t x = 0;
        std::int64_t y = 0;
        std::int64_t g = extended_gcd(absB, absA, x, y);
        if (g != 1) {
            std::cerr << "Unexpected gcd value.\n";
            return 1;
        }
        std::int64_t u = x * (B < 0 ? -1 : 1);
        std::int64_t v = y * (A < 0 ? -1 : 1);

        std::int64_t C = (11 * v + 719 * u) % D;
        if (C < 0) {
            C += D;
        }

        std::int64_t g1 = std::gcd<std::int64_t>(kCycleSeconds, D);
        std::int64_t num_base = kCycleSeconds / g1;
        std::int64_t den_base = D / g1;

        for (std::int64_t n = 0; n < D; ++n) {
            std::int64_t n2 = static_cast<std::int64_t>((__int128)n * C % D);
            if (n2 == n) {
                continue;
            }
            std::int64_t g2 = std::gcd<std::int64_t>(n, den_base);
            std::int64_t num = num_base * (n / g2);
            std::int64_t den = den_base / g2;
            moments.emplace_back(num, den);
        }
    } while (std::next_permutation(idx.begin(), idx.end()));

    std::sort(moments.begin(), moments.end());
    moments.erase(std::unique(moments.begin(), moments.end()), moments.end());

    std::cout << moments.size() << '\n';
    return 0;
}
