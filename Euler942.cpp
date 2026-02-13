#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

constexpr std::uint64_t kOutMod = 1'000'000'007ULL;

std::uint64_t solve(std::uint32_t q) {
    const std::uint32_t half = (q - 1U) / 2U;

    std::vector<std::uint8_t> is_residue(static_cast<std::size_t>(q), 0);
    std::uint64_t r = 1;
    for (std::uint32_t x = 1; x <= half; ++x) {
        is_residue[static_cast<std::size_t>(r)] = 1;
        r += static_cast<std::uint64_t>(2U) * x + 1U;
        if (r >= q) {
            r -= q;
        }
    }

    std::uint64_t pow2 = 1;
    std::uint64_t sum_nonres = 0;
    std::uint64_t signed_sum = 0;

    for (std::uint32_t n = 1; n < q; ++n) {
        pow2 = (pow2 * 2ULL) % kOutMod;
        if (is_residue[n]) {
            signed_sum += pow2;
            if (signed_sum >= kOutMod) {
                signed_sum -= kOutMod;
            }
        } else {
            if (signed_sum >= pow2) {
                signed_sum -= pow2;
            } else {
                signed_sum += kOutMod - pow2;
            }
            sum_nonres += pow2;
            if (sum_nonres >= kOutMod) {
                sum_nonres -= kOutMod;
            }
        }
    }

    if ((q & 7U) == 1U) {
        return (1ULL + 2ULL * sum_nonres) % kOutMod;
    }
    return signed_sum;
}

void run_validations() {
    assert(solve(5) == 6ULL);
    assert(solve(17) == 47'569ULL);
}

}  // namespace

int main() {
    run_validations();
    constexpr std::uint32_t kQ = 74'207'281U;
    std::cout << solve(kQ) << '\n';
    return 0;
}
