#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

struct Options {
    bool run_checkpoints = true;
};

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

bool is_pandigital_1_to_9(std::uint32_t value) {
    int mask = 0;
    int count = 0;
    while (value > 0) {
        const int d = static_cast<int>(value % 10U);
        if (d == 0) {
            return false;
        }
        const int bit = 1 << d;
        if ((mask & bit) != 0) {
            return false;
        }
        mask |= bit;
        value /= 10U;
        ++count;
    }
    return count == 9 && mask == 0x3FE;
}

std::uint32_t leading_nine_digits(const int n) {
    static const long double log10_phi = std::log10((1.0L + std::sqrt(5.0L)) / 2.0L);
    static const long double log10_sqrt5 = std::log10(std::sqrt(5.0L));

    const long double x = static_cast<long double>(n) * log10_phi - log10_sqrt5;
    const long double fractional = x - std::floor(x);
    long double v = std::pow(10.0L, fractional + 8.0L);

    std::uint32_t lead = static_cast<std::uint32_t>(v + 1e-10L);
    while (lead >= 1000000000U) {
        lead /= 10U;
    }
    return lead;
}

std::uint32_t trailing_nine_digits_at_index(const int target_index) {
    if (target_index <= 2) {
        return 1;
    }

    std::uint32_t a = 1;
    std::uint32_t b = 1;
    constexpr std::uint32_t mod = 1000000000U;

    for (int n = 3; n <= target_index; ++n) {
        const std::uint32_t c = static_cast<std::uint32_t>((static_cast<std::uint64_t>(a) + b) % mod);
        a = b;
        b = c;
    }
    return b;
}

int solve() {
    constexpr std::uint32_t mod = 1000000000U;
    std::uint32_t a = 1;
    std::uint32_t b = 1;

    for (int n = 3;; ++n) {
        const std::uint32_t c = static_cast<std::uint32_t>((static_cast<std::uint64_t>(a) + b) % mod);
        a = b;
        b = c;

        if (!is_pandigital_1_to_9(b)) {
            continue;
        }

        const std::uint32_t lead = leading_nine_digits(n);
        if (is_pandigital_1_to_9(lead)) {
            return n;
        }
    }
}

bool run_checkpoints() {
    if (!is_pandigital_1_to_9(trailing_nine_digits_at_index(541))) {
        std::cerr << "Checkpoint failed for trailing pandigital index 541" << '\n';
        return false;
    }
    if (!is_pandigital_1_to_9(leading_nine_digits(2749))) {
        std::cerr << "Checkpoint failed for leading pandigital index 2749" << '\n';
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

    std::cout << solve() << '\n';
    return 0;
}
