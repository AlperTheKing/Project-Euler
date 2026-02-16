#include <algorithm>
#include <complex>
#include <cstdint>
#include <iostream>

namespace {

using i64 = long long;
using Complex = std::complex<i64>;

constexpr int kLimit = 150;
constexpr int kPrimeCount = 16;
constexpr int kPrimes[kPrimeCount] = {
    5, 13, 17, 29, 37, 41, 53, 61,
    73, 89, 97, 101, 109, 113, 137, 149
};

Complex reps[kLimit];

i64 dfs(int idx, const Complex& value) {
    if (idx == kPrimeCount) {
        const i64 a = std::llabs(value.real());
        const i64 b = std::llabs(value.imag());
        return std::min(a, b);
    }
    return dfs(idx + 1, value) +
           dfs(idx + 1, value * reps[kPrimes[idx]]) +
           dfs(idx + 1, value * std::conj(reps[kPrimes[idx]]));
}

void build_representations() {
    for (int i = 1; i * i < kLimit; i += 2) {
        for (int j = 2; i * i + j * j < kLimit; j += 2) {
            reps[i * i + j * j] = Complex(i, j);
        }
    }
}

i64 solve() {
    build_representations();
    i64 sum = 0;
    for (int k = 0; k < kPrimeCount; ++k) {
        sum += dfs(k + 1, reps[kPrimes[k]]);
    }
    return sum;
}

}  // namespace

int main() {
    std::cout << solve() << '\n';
    return 0;
}
