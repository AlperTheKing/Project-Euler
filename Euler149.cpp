#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Options {
    int size = 2000;
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
        if (parse_int_after_prefix(arg, "--size=", options.size)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.size >= 1;
}

std::vector<int> generate_sequence(const int total) {
    std::vector<int> s(static_cast<std::size_t>(total + 1), 0);

    for (int k = 1; k <= std::min(55, total); ++k) {
        const long long kk = static_cast<long long>(k);
        const long long value = (100003LL - 200003LL * kk + 300007LL * kk * kk * kk) % 1000000LL;
        s[static_cast<std::size_t>(k)] = static_cast<int>(value) - 500000;
    }

    for (int k = 56; k <= total; ++k) {
        const int value = (s[static_cast<std::size_t>(k - 24)] + s[static_cast<std::size_t>(k - 55)] + 1000000) % 1000000;
        s[static_cast<std::size_t>(k)] = value - 500000;
    }

    return s;
}

int kadane_extend(const int value, const int current) {
    return std::max(value, current + value);
}

int solve(const int size) {
    const int total = size * size;
    const std::vector<int> s = generate_sequence(total);

    int best = std::numeric_limits<int>::min();

    auto at = [&](const int r, const int c) -> int {
        return s[static_cast<std::size_t>(r * size + c + 1)];
    };

    for (int r = 0; r < size; ++r) {
        int current = std::numeric_limits<int>::min() / 4;
        for (int c = 0; c < size; ++c) {
            current = kadane_extend(at(r, c), current);
            best = std::max(best, current);
        }
    }

    for (int c = 0; c < size; ++c) {
        int current = std::numeric_limits<int>::min() / 4;
        for (int r = 0; r < size; ++r) {
            current = kadane_extend(at(r, c), current);
            best = std::max(best, current);
        }
    }

    for (int start = 0; start < size; ++start) {
        int current = std::numeric_limits<int>::min() / 4;
        for (int r = 0, c = start; r < size && c < size; ++r, ++c) {
            current = kadane_extend(at(r, c), current);
            best = std::max(best, current);
        }
    }
    for (int start = 1; start < size; ++start) {
        int current = std::numeric_limits<int>::min() / 4;
        for (int r = start, c = 0; r < size && c < size; ++r, ++c) {
            current = kadane_extend(at(r, c), current);
            best = std::max(best, current);
        }
    }

    for (int start = 0; start < size; ++start) {
        int current = std::numeric_limits<int>::min() / 4;
        for (int r = 0, c = start; r < size && c >= 0; ++r, --c) {
            current = kadane_extend(at(r, c), current);
            best = std::max(best, current);
        }
    }
    for (int start = 1; start < size; ++start) {
        int current = std::numeric_limits<int>::min() / 4;
        for (int r = start, c = size - 1; r < size && c >= 0; ++r, --c) {
            current = kadane_extend(at(r, c), current);
            best = std::max(best, current);
        }
    }

    return best;
}

bool run_checkpoints() {
    const std::vector<int> s = generate_sequence(100);
    if (s[10] != -393027) {
        std::cerr << "Checkpoint failed for s10" << '\n';
        return false;
    }
    if (s[100] != 86613) {
        std::cerr << "Checkpoint failed for s100" << '\n';
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

    std::cout << solve(options.size) << '\n';
    return 0;
}
