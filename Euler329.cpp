#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using boost::multiprecision::cpp_int;

struct Options {
    int squares = 500;
    std::string sequence = "PPPPNNPPPNPPNPN";
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_string_after_prefix(const std::string& arg, const std::string& prefix, std::string& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    value = arg.substr(prefix.size());
    return !value.empty();
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--squares=", options.squares) ||
            parse_string_after_prefix(arg, "--sequence=", options.sequence)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.squares < 2 || options.sequence.empty()) {
        return false;
    }
    for (char c : options.sequence) {
        if (c != 'P' && c != 'N') {
            return false;
        }
    }
    return true;
}

std::vector<std::uint8_t> sieve_prime_flags(const int limit) {
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit + 1), 1U);
    is_prime[0] = 0U;
    is_prime[1] = 0U;
    for (int p = 2; static_cast<long long>(p) * p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] == 0U) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            is_prime[static_cast<std::size_t>(q)] = 0U;
        }
    }
    return is_prime;
}

cpp_int gcd_cpp(cpp_int a, cpp_int b) {
    while (b != 0) {
        cpp_int t = a % b;
        a = b;
        b = t;
    }
    return a;
}

std::pair<cpp_int, cpp_int> probability_fraction(const int squares, const std::string& sequence) {
    const auto is_prime = sieve_prime_flags(squares);
    std::vector<cpp_int> state(static_cast<std::size_t>(squares + 1), 1);
    cpp_int denom = squares;

    for (std::size_t step = 0; step < sequence.size(); ++step) {
        const char croak = sequence[step];
        for (int pos = 1; pos <= squares; ++pos) {
            const bool prime_square = (is_prime[static_cast<std::size_t>(pos)] != 0U);
            const bool good = (croak == 'P') ? prime_square : !prime_square;
            if (good) {
                state[static_cast<std::size_t>(pos)] *= 2;
            }
        }
        denom *= 3;

        if (step + 1 < sequence.size()) {
            std::vector<cpp_int> next(static_cast<std::size_t>(squares + 1), 0);
            for (int pos = 1; pos <= squares; ++pos) {
                const cpp_int& w = state[static_cast<std::size_t>(pos)];
                if (pos == 1) {
                    next[2] += w * 2;
                } else if (pos == squares) {
                    next[static_cast<std::size_t>(squares - 1)] += w * 2;
                } else {
                    next[static_cast<std::size_t>(pos - 1)] += w;
                    next[static_cast<std::size_t>(pos + 1)] += w;
                }
            }
            denom *= 2;
            state.swap(next);
        }
    }

    cpp_int numer = 0;
    for (int pos = 1; pos <= squares; ++pos) {
        numer += state[static_cast<std::size_t>(pos)];
    }
    const cpp_int g = gcd_cpp(numer, denom);
    numer /= g;
    denom /= g;
    return {numer, denom};
}

bool run_checkpoints() {
    {
        const auto p = probability_fraction(500, "P");
        if (p.first != 119 || p.second != 300) {
            std::cerr << "Checkpoint failed for one-croak probability P" << '\n';
            return false;
        }
    }
    {
        const auto p = probability_fraction(500, "P");
        const auto n = probability_fraction(500, "N");
        if (p.first * n.second + n.first * p.second != p.second * n.second) {
            std::cerr << "Checkpoint failed for P/N complementarity" << '\n';
            return false;
        }
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
    const auto [numer, denom] = probability_fraction(options.squares, options.sequence);
    std::cout << numer << "/" << denom << '\n';
    return 0;
}
