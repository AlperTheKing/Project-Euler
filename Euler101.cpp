#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

namespace {

using i128 = __int128_t;

struct Options {
    int degree = 10;
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

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--degree=", options.degree)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.degree >= 1;
}

i128 term_from_coeffs(const i128 n, const std::vector<i128>& coeffs) {
    i128 value = 0;
    i128 power = 1;
    for (std::size_t k = 0; k < coeffs.size(); ++k) {
        value += coeffs[k] * power;
        power *= n;
    }
    return value;
}

i128 binom_i128(i128 n, int k) {
    if (k < 0) {
        return 0;
    }
    if (k == 0) {
        return 1;
    }
    i128 value = 1;
    for (int i = 1; i <= k; ++i) {
        value = value * (n - i + 1) / i;
    }
    return value;
}

i128 sum_of_fits_from_coeffs(const std::vector<i128>& coeffs) {
    const int degree = static_cast<int>(coeffs.size()) - 1;
    std::vector<i128> values;
    values.reserve(static_cast<std::size_t>(degree + 1));
    for (int n = 1; n <= degree + 1; ++n) {
        values.push_back(term_from_coeffs(n, coeffs));
    }

    std::vector<i128> first_diff_entries;
    first_diff_entries.reserve(values.size());

    std::vector<i128> current = values;
    while (!current.empty()) {
        first_diff_entries.push_back(current[0]);
        if (current.size() == 1) {
            break;
        }
        std::vector<i128> next;
        next.reserve(current.size() - 1);
        for (std::size_t i = 1; i < current.size(); ++i) {
            next.push_back(current[i] - current[i - 1]);
        }
        current = std::move(next);
    }

    i128 answer = 0;
    for (int k = 1; k <= degree; ++k) {
        const int n = k + 1;
        i128 predicted = 0;
        for (int j = 0; j < k; ++j) {
            predicted += binom_i128(n - 1, j) * first_diff_entries[static_cast<std::size_t>(j)];
        }

        const i128 actual = values[static_cast<std::size_t>(k)];
        if (predicted != actual) {
            answer += predicted;
        }
    }

    return answer;
}

i128 sum_of_fits_alternating(const int degree) {
    std::vector<i128> coeffs(static_cast<std::size_t>(degree + 1), 0);
    for (int k = 0; k <= degree; ++k) {
        coeffs[static_cast<std::size_t>(k)] = ((k & 1) == 0) ? 1 : -1;
    }
    return sum_of_fits_from_coeffs(coeffs);
}

std::string i128_to_string(i128 v) {
    if (v == 0) {
        return "0";
    }

    bool negative = (v < 0);
    if (negative) {
        v = -v;
    }

    std::string s;
    while (v > 0) {
        const int digit = static_cast<int>(v % 10);
        s.push_back(static_cast<char>('0' + digit));
        v /= 10;
    }
    if (negative) {
        s.push_back('-');
    }
    std::reverse(s.begin(), s.end());
    return s;
}

bool run_checkpoints() {
    const std::vector<i128> cubic = {0, 0, 0, 1};
    if (sum_of_fits_from_coeffs(cubic) != 74) {
        std::cerr << "Checkpoint failed for cubic sequence" << '\n';
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

    std::cout << i128_to_string(sum_of_fits_alternating(options.degree)) << '\n';
    return 0;
}
