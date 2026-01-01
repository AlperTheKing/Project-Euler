#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr std::int64_t kMod = 1234567891LL;

inline std::int64_t mul_mod(std::int64_t a, std::int64_t b) {
    return static_cast<std::int64_t>((__int128)a * b % kMod);
}

std::vector<int> build_inverses(int limit) {
    std::vector<int> inv(limit + 1, 0);
    if (limit >= 1) {
        inv[1] = 1;
    }
    for (int i = 2; i <= limit; ++i) {
        inv[i] = static_cast<int>(kMod - (kMod / i) * 1LL * inv[kMod % i] % kMod);
    }
    return inv;
}

std::int64_t next_d(int n, std::int64_t d_n, std::int64_t d_n1,
                    const std::vector<int>& inv) {
    const std::int64_t n2 = (static_cast<std::int64_t>(n) * n) % kMod;
    const std::int64_t n1 = n + 1;

    std::int64_t term1 = mul_mod(n2, n1);
    term1 = mul_mod(term1, 2LL * n + 3);
    term1 = mul_mod(term1, 16);
    term1 = mul_mod(term1, d_n);

    std::int64_t quad = (2LL * n2 + 4LL * n + 3) % kMod;
    std::int64_t term2 = mul_mod(quad, d_n1);
    term2 = mul_mod(term2, n1);
    term2 = mul_mod(term2, 4);

    std::int64_t numerator = (term1 + term2) % kMod;

    std::int64_t denom_inv = mul_mod(inv[n], inv[n + 2]);
    denom_inv = mul_mod(denom_inv, inv[n + 3]);
    denom_inv = mul_mod(denom_inv, inv[2 * n + 1]);

    return mul_mod(numerator, denom_inv);
}

void validate(std::int64_t d3, std::int64_t d100) {
    if (d3 != 4) {
        std::cerr << "Validation failed for D(3).\n";
        std::exit(1);
    }
    if (d100 != 1172122931LL) {
        std::cerr << "Validation failed for D(100) modulo " << kMod << ".\n";
        std::exit(1);
    }
}

}  // namespace

int main(int argc, char** argv) {
    int target = 10'000'000;
    bool do_validate = true;
    bool target_set = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--no-validate") {
            do_validate = false;
            continue;
        }
        char* end = nullptr;
        long long val = std::strtoll(arg.c_str(), &end, 10);
        if (end && *end == '\0' && !target_set) {
            target = static_cast<int>(std::max(0LL, val));
            target_set = true;
        }
    }

    if (target <= 0) {
        std::cout << 0 << '\n';
        return 0;
    }

    int compute_n = target;
    if (do_validate) {
        compute_n = std::max(compute_n, 100);
    }

    const int inv_limit = 2 * compute_n + 3;
    std::vector<int> inv = build_inverses(inv_limit);

    std::int64_t d_prev = 0;  // D(1)
    std::int64_t d_curr = 2;  // D(2)
    std::int64_t total = 0;
    if (target >= 1) total = d_prev;
    if (target >= 2) total = (total + d_curr) % kMod;

    std::int64_t d3 = -1;
    std::int64_t d100 = -1;

    for (int n = 1; n + 2 <= compute_n; ++n) {
        std::int64_t d_next = next_d(n, d_prev, d_curr, inv);
        const int idx = n + 2;
        if (idx == 3) d3 = d_next;
        if (idx == 100) d100 = d_next;
        if (idx <= target) {
            total += d_next;
            if (total >= kMod) total -= kMod;
        }
        d_prev = d_curr;
        d_curr = d_next;
    }

    if (do_validate) {
        validate(d3, d100);
    }

    std::cout << total % kMod << '\n';
    return 0;
}
