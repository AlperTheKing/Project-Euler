#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

u128 T(const u128 n) {
    return n * (n + 1);
}

u128 H(const u128 n) {
    const u128 m = n * n * (n + 1);
    return m * (m + 1);
}

u128 U(const u128 n) {
    const u128 h = H(n);
    return h * (h + 1);
}

u64 to_u64(const u128 x) {
    return static_cast<u64>(x);
}

void validate() {
    assert(to_u64(T(1)) == 2);
    assert(to_u64(H(1)) == 6);
    assert(to_u64(H(2)) == 156);
    assert(to_u64(U(1)) == 42);
}

}  // namespace

int main() {
    validate();

    constexpr u64 kMod = 1'000'000'000ULL;

    const u128 n = U(1);
    const u128 value = U(n);
    const u64 answer = static_cast<u64>(value % kMod);

    std::cout << answer << '\n';
    return 0;
}
