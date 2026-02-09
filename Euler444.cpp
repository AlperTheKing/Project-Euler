#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

using u64 = std::uint64_t;
using Big = boost::multiprecision::cpp_dec_float_100;

struct Options {
    u64 n = 100'000'000'000'000ULL;
    u64 k = 20ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = static_cast<u64>(std::stoull(tail));
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--k=", options.k)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1ULL;
}

Big harmonic_small(const u64 n) {
    Big h = 0;
    for (u64 i = 1; i <= n; ++i) {
        h += Big(1) / Big(i);
    }
    return h;
}

Big harmonic_large(const u64 n) {
    if (n <= 1'000'000ULL) {
        return harmonic_small(n);
    }

    const Big x = Big(n);
    const Big inv = Big(1) / x;
    const Big inv2 = inv * inv;
    const Big inv4 = inv2 * inv2;
    const Big inv6 = inv4 * inv2;
    const Big inv8 = inv4 * inv4;
    const Big inv10 = inv8 * inv2;
    const Big inv12 = inv10 * inv2;
    const Big inv14 = inv12 * inv2;
    const Big inv16 = inv14 * inv2;
    const Big inv18 = inv16 * inv2;
    const Big inv20 = inv18 * inv2;

    return log(x) + boost::math::constants::euler<Big>() + inv / 2 - inv2 / 12 + inv4 / 120 -
           inv6 / 252 + inv8 / 240 - Big(5) * inv10 / 660 + Big(691) * inv12 / 32760 -
           inv14 / 12 + Big(3617) * inv16 / 8160 - Big(43867) * inv18 / 14364 +
           Big(174611) * inv20 / 6600;
}

Big binomial_n_plus_k_choose_k(const u64 n, const u64 k) {
    Big c = 1;
    for (u64 i = 1; i <= k; ++i) {
        c *= Big(n + i);
        c /= Big(i);
    }
    return c;
}

Big compute_s_k(const u64 n, const u64 k) {
    const Big comb = binomial_n_plus_k_choose_k(n, k);
    const Big h_big = harmonic_large(n + k);
    const Big h_k = harmonic_small(k);
    return comb * (h_big - h_k);
}

std::string normalize_scientific(const std::string& s) {
    const std::size_t pos = s.find('e');
    if (pos == std::string::npos) {
        return s;
    }

    const std::string mantissa = s.substr(0, pos);
    std::string exponent = s.substr(pos + 1);

    bool negative = false;
    if (!exponent.empty() && (exponent[0] == '+' || exponent[0] == '-')) {
        negative = (exponent[0] == '-');
        exponent.erase(exponent.begin());
    }

    while (exponent.size() > 1U && exponent[0] == '0') {
        exponent.erase(exponent.begin());
    }

    return mantissa + (negative ? "e-" : "e") + exponent;
}

std::string format_scientific_10sig(const Big& value) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(9) << value;
    return normalize_scientific(oss.str());
}

bool run_checkpoints() {
    std::ostringstream e111;
    e111 << std::fixed << std::setprecision(4) << harmonic_small(111ULL);
    if (e111.str() != "5.2912") {
        std::cerr << "Checkpoint failed: E(111)\n";
        return false;
    }

    const std::string s3 = format_scientific_10sig(compute_s_k(100ULL, 3ULL));
    if (s3 != "5.983679014e5") {
        std::cerr << "Checkpoint failed: S3(100)\n";
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

    std::cout << format_scientific_10sig(compute_s_k(options.n, options.k)) << '\n';
    return 0;
}
